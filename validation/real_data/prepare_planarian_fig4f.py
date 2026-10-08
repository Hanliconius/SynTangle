"""Import the authors' Fig. 4f BUSCO tables; no BUSCO or solver execution."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
from prepare_published_inputs import download, digest

COMMIT = '91ae879d1d72386d393f38b5e5b4d111a2f1458b'
SPECIES = {
 'schMedS3_h1': 'schMedS3_h1', 'schMedS3_h2': 'schMedS3_h2',
 'schPol2': 'schPol2', 'schNov1': 'schNov1', 'schLug1': 'schLug1',
 'cloSin': 'clonorchis_sinensis.PRJNA386618',
 'schMan': 'schistosoma_mansoni.PRJEA36577',
 'hymMic': 'hymenolepis_microstoma.PRJEB124',
 'taeMul': 'taenia_multiceps.PRJNA307624',
}
def canonical(sp, sequence):
    sequence = re.sub(r':[0-9]+-[0-9]+', '', sequence.replace('manual_scaffold_', 'chr'))
    if sp.startswith('schMedS3_'): sequence = re.sub(r'chr([0-9])_h.', r'chr\1', sequence)
    if sp == 'cloSin':
        if sequence not in list(map(str, range(1, 8))): return None
        sequence = 'chr' + sequence
    if sp == 'hymMic': sequence = re.sub(r'HMN_0(.)_pilon', r'chr\1', sequence)
    if sp == 'schMan':
        if 'H0' in sequence or 'U0' in sequence: return None
        sequence = sequence.replace('SM_V7_', 'chr')
    if sp == 'taeMul': sequence = sequence.replace('LG', 'chr')
    if sp == 'schPol2' and not re.search(r'chr[1234]$', sequence): return None
    if not re.fullmatch(r'chr(?:[1-9][0-9]?|ZW)', sequence): return None
    return sequence
def rank(chrom):
    return 100 if chrom == 'chrZW' else int(chrom[3:])

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--local-author-repo', type=Path)
    a = p.parse_args()
    raw = a.output / 'raw'; raw.mkdir(parents=True, exist_ok=True)
    rows = defaultdict(list); lengths = defaultdict(dict)
    excluded = Counter(); sources = []
    base = f'https://raw.githubusercontent.com/Jeremias-Brand/PlanarianGenomeAnalysis/{COMMIT}/busco_synteny/'
    for sp, folder in SPECIES.items():
        for rel in [f'chromsize/{folder}.chromsize',
                    f'01_busco/{folder}/busco/run_metazoa_odb10/full_table.tsv']:
            dest = raw / sp / Path(rel).name
            dest.parent.mkdir(parents=True, exist_ok=True)
            if a.local_author_repo:
                import shutil
                shutil.copyfile(a.local_author_repo / 'busco_synteny' / rel, dest)
            else: download(base + rel, dest)
            sources.append({'url': base + rel, 'sha256': digest(dest), 'species': sp})
        native_lengths = {}
        for line in (raw / sp / f'{folder}.chromsize').read_text().splitlines():
            seq, size = line.split()[:2]
            native_lengths[seq] = int(size)
        for line in (raw / sp / 'full_table.tsv').read_text().splitlines():
            if not line or line.startswith('#'): continue
            t = line.split('\t')
            if len(t) < 8 or t[1] != 'Complete':
                excluded[sp + ':not_Complete'] += 1; continue
            h, status, native, start, end, strand = t[:6]
            chrom = canonical(sp, native)
            if chrom is None:
                excluded[sp + ':outside_author_chromosomes'] += 1; continue
            length = native_lengths.get(native) or native_lengths.get(re.sub(r':[0-9]+-[0-9]+', '', native.replace('manual_scaffold_', 'chr')))
            if length is None:
                # Authors remove the coordinate suffix when parsing the BUSCO sequence.
                length = native_lengths.get(re.sub(r':[0-9]+-[0-9]+', '', native))
            if length is None: raise ValueError(f'Missing chromosome length {sp}/{native}')
            start, end = sorted((int(start), int(end)))
            if not 0 <= start < end <= length: raise ValueError(f'Invalid source coordinates {sp}/{h}')
            if chrom in lengths[sp] and lengths[sp][chrom] != length:
                raise ValueError(f'Conflicting chromosome lengths {sp}/{chrom}')
            lengths[sp][chrom] = length
            rows[h].append({'sp': sp, 'chr': chrom, 'start': start, 'end': end,
                            'strand': strand, 'length': length, 'native_sequence': native})
    retained = {}
    for h, occurrences in rows.items():
        if len(occurrences) == 9 and {r['sp'] for r in occurrences} == set(SPECIES):
            retained[h] = occurrences
        else: excluded['not_one_Complete_occurrence_in_all_nine'] += 1
    if not retained: raise ValueError('No shared single-copy BUSCOs')
    blocks = defaultdict(list); colours = {}
    for h, occurrences in retained.items():
        colours[h] = next(r['chr'] for r in occurrences if r['sp'] == 'schMan')
        for r in occurrences:
            # Author Fig. 4f uses rounded gene START / chromosome length,
            # equal-width chromosomes and native orientations.
            x = round(r['start'] / r['length'], 3)
            blocks[r['sp'], r['chr']].append({
                'occurrence_id': r['sp'] + ':' + h, 'homology_id': h,
                'start': max(0, x - .000001), 'end': min(1, x + .000001),
                'strand': r['strand']})
    fixture = {'fixture_version': 1, 'id': 'flatworms_fig4f_busco',
        'title': 'Published flatworm Fig. 4f BUSCO positions',
        'purpose': 'Paper BUSCO display reconstruction; no new homology inference',
        'orientation_constraints': [], 'species': [
            {'id': sp, 'chromosomes': [
                {'id': chrom, 'length': 1, 'display_rank': rank(chrom),
                 'blocks': blocks[sp, chrom]} for chrom in sorted(lengths[sp], key=rank)]}
            for sp in SPECIES]}
    from syntangle.fixtures import fixture_from_dict
    fixture_from_dict(fixture)
    links = len(retained) * 8
    provenance = {'paper': 'https://doi.org/10.1038/s41467-024-52380-9',
        'figure': '4f', 'source_git_commit': COMMIT, 'sources': sources,
        'source_kind': 'Authors deposited single-copy BUSCO coordinates, Fig. 4f',
        'retained_links': links, 'retained_buscos': len(retained),
        'pair_audit': [{'species': [x, y], 'links': len(retained)}
                       for x, y in zip(SPECIES, list(SPECIES)[1:])],
        'chromosome_extent': 'Equal-width chromosomes; rounded gene-start fractions, as in author Fig. 4f script',
        'input_layout_note': 'Fig. 4f order: numbered chromosomes then ZW; native orientations; author normalised gene starts.',
        'colour_reference': 'schMan', 'homology_colour_reference': colours,
        'scope': 'Reconstruction from published BUSCO tables and author script, not panel-b GENESPACE blocks.',
        'exclusions': dict(excluded), 'original_coordinates': retained}
    a.output.mkdir(parents=True, exist_ok=True)
    (a.output / 'fixture.json').write_text(json.dumps(fixture, indent=2) + '\n')
    (a.output / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(f'PREPARED FIG4F: {len(retained)} BUSCOs; {links} adjacent links; ' +
          f'{sum(len(s["chromosomes"]) for s in fixture["species"])} chromosomes', flush=True)
if __name__ == '__main__': main()
