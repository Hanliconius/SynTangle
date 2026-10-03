from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def audit_groups(
    rows: list[dict[str, str]],
    *,
    require_proven: bool,
) -> dict[str, object]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        biology_id = row.get("biology_id", "").strip()
        if biology_id:
            grouped[biology_id].append(row)

    paired = {
        biology_id: group
        for biology_id, group in grouped.items()
        if len(group) > 1
    }
    if not paired:
        raise ValueError("No paired biology_id groups found")

    details = []
    failures = []

    for biology_id, group in sorted(paired.items()):
        modes = tuple(sorted(row["tangle_mode"] for row in group))
        fingerprints = {
            row["biological_fingerprint"]
            for row in group
        }
        native = {
            row.get("hidden_native_crossings", "")
            for row in group
        }
        optima = {
            row["optimized_crossings"]
            for row in group
            if row["optimality_status"] == "proven optimum"
        }
        statuses = {
            row["optimality_status"]
            for row in group
        }
        all_proven = statuses == {"proven optimum"}
        same_biology = len(fingerprints) == 1
        same_native = len(native) == 1
        same_optimum = all_proven and len(optima) == 1

        if not same_biology:
            failures.append(
                f"{biology_id}: biological fingerprints differ across tangle modes"
            )
        if not same_native:
            failures.append(
                f"{biology_id}: hidden native crossing baselines differ"
            )
        if require_proven and not all_proven:
            failures.append(
                f"{biology_id}: not every tangle presentation was proven optimal"
            )
        if all_proven and not same_optimum:
            failures.append(
                f"{biology_id}: proven optimum differs across presentations"
            )

        details.append(
            {
                "biology_id": biology_id,
                "modes": list(modes),
                "same_biology": same_biology,
                "same_native": same_native,
                "all_proven": all_proven,
                "same_proven_optimum": same_optimum,
                "optimized_crossings": sorted(
                    {
                        int(row["optimized_crossings"])
                        for row in group
                    }
                ),
                "initial_crossings": {
                    row["tangle_mode"]: int(row["initial_crossings"])
                    for row in group
                },
                "wall_seconds": {
                    row["tangle_mode"]: float(row["wall_seconds"])
                    for row in group
                },
            }
        )

    return {
        "paired_group_count": len(details),
        "all_same_biology": all(
            item["same_biology"] for item in details
        ),
        "all_same_native": all(
            item["same_native"] for item in details
        ),
        "all_proven": all(
            item["all_proven"] for item in details
        ),
        "all_same_proven_optimum": all(
            item["same_proven_optimum"] for item in details
        ),
        "failures": failures,
        "groups": details,
    }


def write_markdown(audit: dict[str, object], path: Path) -> None:
    lines = [
        "# Paired presentation-tangle invariance",
        "",
        f"- Paired biological groups: {audit['paired_group_count']}",
        f"- Identical biological evidence within every group: {audit['all_same_biology']}",
        f"- Identical hidden native baseline within every group: {audit['all_same_native']}",
        f"- Every presentation proven optimal: {audit['all_proven']}",
        f"- Proven optimum invariant to presentation: {audit['all_same_proven_optimum']}",
        "",
        "| Biology | Initial C by tangle | Proven C* | Status |",
        "|---|---|---:|---|",
    ]

    for item in audit["groups"]:
        initial = ", ".join(
            f"{mode}={value}"
            for mode, value in sorted(item["initial_crossings"].items())
        )
        optima = ", ".join(
            str(value) for value in item["optimized_crossings"]
        )
        status = (
            "PASS"
            if (
                item["same_biology"]
                and item["same_native"]
                and item["all_proven"]
                and item["same_proven_optimum"]
            )
            else "CHECK"
        )
        lines.append(
            f"| {item['biology_id']} | {initial} | {optima} | {status} |"
        )

    if audit["failures"]:
        lines.extend(["", "## Failures", ""])
        lines.extend(f"- {failure}" for failure in audit["failures"])

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("results_tsv")
    parser.add_argument("--json")
    parser.add_argument("--markdown")
    parser.add_argument(
        "--require-proven",
        action="store_true",
        help="Fail if any presentation in a paired group lacks a proof.",
    )
    args = parser.parse_args()

    audit = audit_groups(
        read_tsv(Path(args.results_tsv)),
        require_proven=args.require_proven,
    )

    if args.json:
        Path(args.json).write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if args.markdown:
        write_markdown(audit, Path(args.markdown))

    print(
        "paired_groups={paired_group_count} same_biology={all_same_biology} "
        "all_proven={all_proven} invariant_optimum={all_same_proven_optimum}".format(
            **audit
        )
    )

    if audit["failures"]:
        for failure in audit["failures"]:
            print("FAIL:", failure)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
