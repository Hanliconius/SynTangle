"""Public-evidence layout comparison; no ancestral/hidden files are read."""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import random
import subprocess
import time
from pathlib import Path

from syntangle import LayoutState, load_validation_bundle, optimize_auto, score_crossings
from syntangle.layout import initial_layout_state
from syntangle.visualize import render_layout_state_svg


def read_tsv(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def validate_state(fixture, state):
    expected = set(fixture.chromosome_refs)
    flattened = [ref for refs in state.chromosome_order.values() for ref in refs]
    if set(flattened) != expected or len(flattened) != len(expected):
        raise ValueError("Layout lost or duplicated chromosomes")
    for species in fixture.species_ids:
        if any(ref.species_id != species for ref in state.chromosome_order[species]):
            raise ValueError("Layout reassigned species")
    if set(state.chromosome_orientation) != expected:
        raise ValueError("Layout orientation is incomplete")
    if any(value not in (-1, 1) for value in state.chromosome_orientation.values()):
        raise ValueError("Illegal orientation")


def improve_flips(fixture, start, seed, restarts=3):
    """Greedy flip assistance; fixed GENESPACE order, no ancestry or exact solver."""
    rng = random.Random(seed)
    refs = sorted(fixture.chromosome_refs)
    best = start
    best_score = score_crossings(fixture, best).crossings
    evaluations = 0
    for restart in range(restarts):
        orientation = dict(start.chromosome_orientation)
        if restart:
            orientation = {ref: rng.choice((-1, 1)) for ref in refs}
        state = LayoutState(start.chromosome_order, orientation)
        score = score_crossings(fixture, state).crossings
        while True:
            chosen = None
            chosen_score = score
            for ref in refs:
                signs = dict(state.chromosome_orientation)
                signs[ref] *= -1
                candidate = LayoutState(start.chromosome_order, signs)
                value = score_crossings(fixture, candidate).crossings
                evaluations += 1
                if value < chosen_score:
                    chosen, chosen_score = candidate, value
            if chosen is None:
                break
            state, score = chosen, chosen_score
        if score < best_score:
            best, best_score = state, score
    return best, evaluations


def export_native_input(fixture, directory):
    directory.mkdir(parents=True, exist_ok=True)
    initial = initial_layout_state(fixture)
    ranks = {ref: rank for refs in initial.chromosome_order.values()
             for rank, ref in enumerate(refs, 1)}
    bed, clens = [], []
    for chromosome in fixture.chromosomes:
        ref = chromosome.ref
        clens.append(dict(genome=ref.species_id, chr=ref.chromosome_id,
                          ordByFun=ranks[ref], length=chromosome.length))
        for block in chromosome.blocks:
            midpoint = (block.start + block.end) / 2
            if initial.chromosome_orientation[ref] == -1:
                midpoint = chromosome.length - midpoint
            bed.append(dict(genome=ref.species_id, chr=ref.chromosome_id,
                            ord=midpoint, og=block.homology_id,
                            noAnchor="FALSE", isArrayRep="TRUE"))
    write_tsv(directory / "bed.tsv", bed)
    write_tsv(directory / "clens.tsv", clens)
    write_tsv(directory / "species.tsv", [dict(species_id=s) for s in fixture.species_ids])


def load_native_states(fixture, path):
    lookup = {(ref.species_id, ref.chromosome_id): ref for ref in fixture.chromosome_refs}
    initial = initial_layout_state(fixture)
    groups = {}
    for row in read_tsv(path):
        groups.setdefault(row["variant"], []).append(row)
    for variant, rows in groups.items():
        order = {}
        for species in fixture.species_ids:
            selected = sorted((row for row in rows if row["genome"] == species),
                              key=lambda row: int(row["plotOrd"]))
            order[species] = tuple(lookup[(species, row["chr"])] for row in selected)
        state = LayoutState(order, dict(initial.chromosome_orientation))
        validate_state(fixture, state)
        yield variant, rows[0], state


def write_case_report(output, fixture, results, states):
    cards = []
    for row in results:
        method = row["method"]
        if method not in states:
            continue
        state = states[method]
        svg = render_layout_state_svg(fixture, state, title=method,
                                      crossing_count=row["crossings"])
        (output / f"{method}.svg").write_text(svg)
        (output / f"{method}.layout.json").write_text(json.dumps(state.to_dict(), indent=2))
        cards.append(f'<section><p>{html.escape(method)}: {row["crossings"]} crossings; '
                     f'{row["seconds"]:.4f} s</p>{svg}</section>')
    (output / "comparison.html").write_text(
        '<!doctype html><meta charset="utf-8"><title>Layout comparison</title>'
        '<style>body{font:16px sans-serif;margin:24px;background:#f5f5f5}'
        '.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}'
        'section{background:white;padding:12px}svg{width:100%;height:auto}'
        '@media(max-width:850px){.grid{grid-template-columns:1fr}}</style>'
        f'<h1>{html.escape(fixture.fixture_id)}</h1>'
        '<p>Same physical coordinates, homologies, species rows and renderer. '
        'Whole chromosomes only. GENESPACE results below are the best over all '
        'species as references and syntenyWeight 1/0.5; timings include all trials. '
        'Flip assistance is three-start greedy search with chromosome orders fixed. '
        'Syntangle may change all chromosome orders; GENESPACE fixes its chosen reference order. '
        'This tests layout methods, not GENESPACE homology discovery or full native plotting.</p>'
        '<div class="grid">' + "".join(cards) + '</div>'
    )


def run_case(args):
    root = Path(args.root).resolve()
    manifest = read_tsv(root / "benchmark_manifest.tsv")
    entry = manifest[args.index]
    case = root / entry["case_dir"]
    output = root / "comparison" / entry["case_id"]
    output.mkdir(parents=True, exist_ok=True)
    for stale in ("COMPLETE", "results.tsv", "solver_result.json", "solver_started.json"):
        (output / stale).unlink(missing_ok=True)
    fixture = load_validation_bundle(case)
    results, states = [], {}
    fingerprint = hashlib.sha256(b"".join((case / name).read_bytes() for name in (
        "species.tsv", "chromosomes.tsv", "occurrences.tsv", "input_display_state.tsv"
    ))).hexdigest()
    def add(method, state, seconds, status="heuristic", **extra):
        validate_state(fixture, state)
        states[method] = state
        results.append(dict(case_id=entry["case_id"], method=method,
                            crossings=score_crossings(fixture, state).crossings,
                            seconds=seconds, status=status,
                            species_count=entry["species_count"],
                            ancestor_chromosomes=entry["ancestor_chromosomes"],
                            events_per_branch=entry["events_per_branch"],
                            fingerprint=fingerprint, details=json.dumps(extra, sort_keys=True)))
        write_tsv(output / "results.tsv", results)
        write_case_report(output, fixture, results, states)
        print(f'{entry["case_id"]} {method}: C={results[-1]["crossings"]}, '
              f'{seconds:.4f}s, {status}', flush=True)

    initial = initial_layout_state(fixture)
    add("input", initial, 0, "input")
    export_native_input(fixture, output / "public_native_input")
    helper = Path(__file__).with_name("genespace_native_order.R")
    process_started = time.perf_counter()
    subprocess.run([args.micromamba, "run", "-n", args.genespace_env, "Rscript",
                    str(helper), str(output / "public_native_input"), str(output)],
                   check=True, timeout=180)
    native_process_seconds = time.perf_counter() - process_started
    variants = list(load_native_states(fixture, output / "native_orders.tsv"))
    if not variants:
        raise ValueError("No native GENESPACE layouts")
    native_seconds = sum(float(row["ordering_seconds"]) for _, row, _ in variants)
    best_native = min(variants, key=lambda v: (score_crossings(fixture, v[2]).crossings, v[0]))
    add("GENESPACE", best_native[2], native_seconds,
        selected_variant=best_native[0], native_process_seconds=native_process_seconds)
    flip_seconds = 0.0
    assisted, variants_audit = [], []
    for number, (variant, row, state) in enumerate(variants):
        flip_started = time.perf_counter()
        flipped, evaluations = improve_flips(fixture, state, int(entry["seed"]) + number)
        flip_seconds += time.perf_counter() - flip_started
        assisted.append((variant, flipped))
        variants_audit.append(dict(variant=variant, reference=row["reference"],
                                   synteny_weight=row["synteny_weight"],
                                   raw_crossings=score_crossings(fixture, state).crossings,
                                   assisted_crossings=score_crossings(fixture, flipped).crossings,
                                   flip_candidate_evaluations=evaluations))
        (output / f"{variant}.layout.json").write_text(json.dumps(state.to_dict(), indent=2))
        (output / f"{variant}.flips.layout.json").write_text(json.dumps(flipped.to_dict(), indent=2))
    best_assisted = min(assisted, key=lambda v: (score_crossings(fixture, v[1]).crossings, v[0]))
    write_tsv(output / "reference_variants.tsv", variants_audit)
    add("GENESPACE_plus_flips", best_assisted[1], native_seconds + flip_seconds,
        selected_variant=best_assisted[0], native_process_seconds=native_process_seconds,
        flip_seconds=flip_seconds, restarts=3)

    # Checkpoint the baselines before entering the potentially expensive solver.
    (output / "solver_started.json").write_text(json.dumps(dict(
        transition_cap=args.transition_cap, branch_node_cap=args.branch_node_cap,
        started_at=time.time(), fingerprint=fingerprint)))
    started = time.perf_counter()
    result = optimize_auto(fixture, transition_cap_per_component=args.transition_cap,
                           branch_node_cap_per_component=args.branch_node_cap,
                           local_restarts=2, component_workers=1, seed=int(entry["seed"]))
    seconds = time.perf_counter() - started
    status = result.layout.optimality_status
    add("Syntangle", result.layout.optimized_state, seconds, status,
        solver=result.solver, solver_details=result.details)
    (output / "solver_result.json").write_text(json.dumps(result.to_dict(), indent=2))
    if "proven" in status and results[-1]["crossings"] > results[-2]["crossings"]:
        raise AssertionError("GENESPACE beat a claimed proven optimum: solver regression")
    (output / "COMPLETE").write_text("PASS\n")


def collect(root):
    root = Path(root).resolve()
    rows, comparisons = [], []
    for entry in read_tsv(root / "benchmark_manifest.tsv"):
        output = root / "comparison" / entry["case_id"]
        path = output / "results.tsv"
        current = read_tsv(path) if path.exists() else []
        rows.extend(current)
        by_method = {r["method"]: r for r in current}
        scores = [by_method.get(m, {}).get("crossings", "missing") for m in (
            "input", "GENESPACE", "GENESPACE_plus_flips", "Syntangle")]
        complete = (output / "COMPLETE").exists()
        outcome = "incomplete; inspect task log"
        if complete:
            gs, st = map(int, (scores[2], scores[3]))
            outcome = "tie" if gs == st else ("Syntangle fewer" if st < gs else "GENESPACE+flips fewer")
        comparisons.append((entry["case_id"], scores, outcome))
    if rows:
        write_tsv(root / "comparison_results.tsv", rows)
    lines = ["# GENESPACE layout pilot", "", "| Case | Input C | GS C | GS + flips C | Syntangle C | Outcome |",
             "|---|---:|---:|---:|---:|---|"]
    for case, scores, outcome in comparisons:
        lines.append("| " + " | ".join([case, *map(str, scores), outcome]) + " |")
    lines.extend(["", "Missing Syntangle results are not wins for either method.",
                  "Only COMPLETE cases are classified. Check Slurm exit/timeout logs for incomplete cases."])
    (root / "comparison_summary.md").write_text("\n".join(lines) + "\n")
    links = []
    for case, scores, outcome in comparisons:
        links.append(f'<tr><td><a href="comparison/{case}/comparison.html">{html.escape(case)}</a></td>'
                     + "".join(f'<td>{value}</td>' for value in scores)
                     + f'<td>{html.escape(outcome)}</td></tr>')
    (root / "comparison_index.html").write_text(
        '<!doctype html><meta charset="utf-8"><title>GENESPACE pilot</title>'
        '<style>body{font:16px sans-serif;margin:24px}td,th{padding:8px;text-align:left}</style>'
        '<h1>GENESPACE layout pilot</h1><p>Click a case for matched vector figures. '
        'This is a small synthetic pilot, not evidence of general superiority.</p>'
        '<table><tr><th>Case</th><th>Input C</th><th>GS C</th><th>GS + flips C</th>'
        '<th>Syntangle C</th><th>Outcome</th></tr>' + "".join(links) + '</table>')
    print("\n".join(lines))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root")
    parser.add_argument("--index", type=int)
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--micromamba", default="micromamba")
    parser.add_argument("--genespace-env", default="lep_busco_painter_clean")
    parser.add_argument("--transition-cap", type=int, default=100000)
    parser.add_argument("--branch-node-cap", type=int, default=25000)
    args = parser.parse_args()
    if args.collect:
        collect(args.root)
    elif args.index is not None:
        run_case(args)
    else:
        parser.error("Specify --index or --collect")


if __name__ == "__main__":
    main()
