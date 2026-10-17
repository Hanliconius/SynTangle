"""Strict import of directed pairwise links without filtering or coordinate edits."""
import csv
import math
from pathlib import Path
from .fixtures import fixture_from_dict

CHROM_COLUMNS = 'species chromosome length rank orientation'.split()
LINK_COLUMNS = 'link_id species1 chromosome1 start1 end1 species2 chromosome2 start2 end2'.split()


def read_table(path, columns):
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not set(columns) <= set(reader.fieldnames or []):
            raise ValueError(f'{path}: required columns: {", ".join(columns)}')
        rows = list(reader)
    if not rows or any(any(r.get(c) in (None, '') for c in columns) for r in rows):
        raise ValueError(f'{path}: empty table or missing required values')
    return rows


def load_tables(chromosomes, links):
    species = {}
    refs = {}
    for r in read_table(chromosomes, CHROM_COLUMNS):
        sp, ch = r['species'], r['chromosome']
        if ':' in sp or ':' in ch:
            raise ValueError('Species and chromosome IDs must not contain colons')
        key = sp, ch
        if key in refs:
            raise ValueError(f'Duplicate chromosome: {key}')
        length = float(r['length'])
        if not math.isfinite(length) or length <= 0:
            raise ValueError(f'Invalid length: {key}')
        row = dict(id=ch, length=length, display_rank=int(r['rank']),
                   display_orientation=int(r['orientation']), blocks=[])
        refs[key] = row
        species.setdefault(sp, []).append(row)
    order = list(species)
    seen = set()
    for i, r in enumerate(read_table(links, LINK_COLUMNS)):
        if r['link_id'] in seen:
            raise ValueError('Duplicate link_id: ' + r['link_id'])
        seen.add(r['link_id'])
        a, b = r['species1'], r['species2']
        if a not in order or b not in order or abs(order.index(a)-order.index(b)) != 1:
            raise ValueError('Links must join adjacent displayed species: ' + r['link_id'])
        for side in (1, 2):
            key = r[f'species{side}'], r[f'chromosome{side}']
            if key not in refs:
                raise ValueError(f'Link references missing chromosome: {key}')
            start, end = float(r[f'start{side}']), float(r[f'end{side}'])
            if not all(map(math.isfinite, (start, end))) or start == end:
                raise ValueError('Endpoints must be distinct and finite: ' + r['link_id'])
            refs[key]['blocks'].append(dict(occurrence_id=f'link{i}:{side}',
                homology_id=r['link_id'], start=min(start,end), end=max(start,end),
                strand='+' if start < end else '-'))
    data = dict(fixture_version=1, id='user_links', title='Supplied synteny links',
                purpose='Optimize fixed pairwise evidence',
                species=[dict(id=sp, chromosomes=rows) for sp,rows in species.items()])
    return fixture_from_dict(data), data
