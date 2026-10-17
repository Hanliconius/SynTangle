"""Capture the exact Bluesky target and audit publicly deposited source inputs.

Run on Pegasus. This does not infer biological links from raster pixels.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import zipfile

from prepare_published_inputs import download, digest

PAPER = 'https://www.nature.com/articles/s41586-026-11057-7'
POST = 'https://bsky.app/profile/hanliconius.bsky.social/post/3mw4xmy2sgs2q'
AUTHOR_COMMIT = 'eae73b27ea4730da6a26477a8edf6dce79f3e011'
DEPOSITS = ['10.17617/3.IXRT5Y', '10.17617/3.DGOJRQ']
ASSETS = {
    333736: ('R_tenuis_REC_h1', 'Rhync_tenuis_REC.hap1.chr.fasta'),
    333735: ('R_tenuis_REC_h1', 'Rhync_tenuis_REC.hap1.chr.helixer.gff3'),
    333744: ('R_austrobrasiliensis_paper_deposit', 'Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta'),
    333741: ('R_austrobrasiliensis_paper_deposit', 'Rhync_austrobrasiliensis_6344D.hic.hap1.chr.helixer.gff3'),
    336260: ('R_breviuscula_companion_h1', 'Rhync_breviuscula.hap1.chr.fasta'),
}


def inventory_sequences(path):
    rows, name, length = [], None, 0
    with path.open() as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    rows.append({'sequence': name, 'length': length})
                name, length = line[1:].split()[0], 0
            else:
                length += len(line.strip())
    if name is not None:
        rows.append({'sequence': name, 'length': length})
    if not rows or any(r['length'] <= 0 for r in rows):
        raise ValueError(f'Invalid FASTA: {path}')
    return rows


def published_geometry(page):
    bars = []
    for d in page.get_drawings():
        r = d['rect']
        if (d['type'] == 'fs' and d['fill'] == (1., 1., 1.)
                and 115 < r.x0 < r.x1 < 296 and 65 < r.y0 < r.y1 < 146):
            bars.append(r)
    rows = []
    for sp, centre, count in [('R_tenuis_REC_h1', 70.7, 2),
                              ('R_austrobrasiliensis', 105.1, 3),
                              ('R_breviuscula', 139.5, 5)]:
        selected = sorted((r for r in bars if abs((r.y0+r.y1)/2-centre) < 1), key=lambda r:r.x0)
        if len(selected) != count:
            raise ValueError(f'Published chromosome geometry changed: {sp}')
        for number, r in enumerate(selected, 1):
            rows.append({'species':sp, 'published_chromosome_label':number,
                         'x0_pdf_points':r.x0, 'x1_pdf_points':r.x1,
                         'y0_pdf_points':r.y0, 'y1_pdf_points':r.y1})
    return rows


def main():
    import fitz
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True, type=Path)
    p.add_argument('--download-assemblies', action='store_true')
    args = p.parse_args()
    out = args.run.resolve(); src = out/'sources'; src.mkdir(parents=True, exist_ok=True)
    download(PAPER+'.pdf', src/'published.pdf')
    download(f'https://raw.githubusercontent.com/Raina-M/Rhynchospora_tenuis_project/{AUTHOR_COMMIT}/run_genespace.R', src/'run_genespace.R')
    download('https://public.api.bsky.app/xrpc/app.bsky.feed.getPostThread?uri=at://did:plc:b3t6tkvujvmh5t7d3dhok52t/app.bsky.feed.post/3mw4xmy2sgs2q&depth=6',src/'bluesky_thread.json')
    doc = fitz.open(src/'published.pdf'); page = doc[2]
    if 'REC haplotype' not in page.get_text():
        raise ValueError('Published figure page changed; inspect before cropping')
    geometry = published_geometry(page)
    with (out/'published_chromosome_geometry.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(geometry[0]),delimiter='\t');w.writeheader();w.writerows(geometry)
    panel = fitz.open(); panel.new_page(width=260,height=108)
    panel[0].show_pdf_page(fitz.Rect(0,0,260,108),doc,2,clip=fitz.Rect(39,49,299,157))
    panel.save(out/'published_fig1a.pdf'); panel.close()
    records, sequences = [], {}
    for doi in DEPOSITS:
        path=src/(doi.split('.')[-1]+'.inventory.json')
        download(f'https://edmond.mpg.de/api/datasets/:persistentId/?persistentId=doi:{doi}',path)
        data=json.loads(path.read_text())
        if data.get('status')!='OK':raise ValueError(f'Deposit unavailable: {doi}')
        for entry in data['data']['latestVersion']['files']:
            d=entry['dataFile']; records.append({'doi':doi,'id':d['id'],'filename':d['filename'],
                'bytes':d['filesize'],'restricted':entry.get('restricted',False),'checksum':d.get('checksum')})
            if args.download_assemblies and d['id'] in ASSETS:
                group, filename=ASSETS[d['id']]
                if filename!=d['filename']:raise ValueError('Deposit file identity changed')
                if entry.get('restricted'):raise ValueError('Selected file now restricted')
                target=src/group/filename
                download(f"https://edmond.mpg.de/api/access/datafile/{d['id']}",target)
                if target.stat().st_size!=d['filesize']:raise ValueError(f'File size mismatch: {target}')
                h=hashlib.md5()
                with target.open('rb') as f:
                    for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
                if h.hexdigest()!=d['checksum']['value']:raise ValueError(f'Deposit checksum mismatch: {target}')
                if filename.endswith('.fasta'):sequences[group]=inventory_sequences(target)
    audit={'paper':PAPER,'panel':'1a','bluesky_post':POST,'author_commit':AUTHOR_COMMIT,
        'published_order':{'R_tenuis_REC_h1':[1,2],'R_austrobrasiliensis':[1,2,3],'R_breviuscula':[1,2,3,4,5]},
        'geometry_units':'PDF points, not inferred base-pair coordinates',
        'ribbons_in_published_pdf':'raster image; original biological endpoints not recoverable exactly from PDF',
        'original_genespace_block_table_found':False,'syntangle_comparison_complete':False,
        'sequence_inventories':sequences,'deposited_files':records,
        'limitations':['Companion breviuscula assembly identity must be checked against paper inputs.',
                       'Austrobrasiliensis is polyploid; do not silently treat all deposited sequences as the three plotted chromosomes.',
                       'Author plotting script covers twenty haplotypes and differs from final three-row published panel.'],
        'source_sha256':{str(x.relative_to(out)):digest(x) for x in src.rglob('*') if x.is_file()}}
    (out/'source_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    (out/'ATTRIBUTION.txt').write_text('Zhang et al., Nature (2026), DOI 10.1038/s41586-026-11057-7. Published Figure 1a cropped without changing its plotted evidence. Article CC BY 4.0. No optimization claim is made by this source capture.\n')
    with zipfile.ZipFile(out/'SynTangle_rhynchospora_sources.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
        for x in sorted(out.rglob('*')):
            if x.is_file() and x.suffix!='.zip' and 'tools' not in x.relative_to(out).parts:
                z.write(x,x.relative_to(out))
    print('EXACT PUBLISHED PANEL:',out/'published_fig1a.pdf',flush=True)
    print('SOURCE AUDIT:',out/'source_audit.json',flush=True)
    print('Original block coordinates still required; no SynTangle comparison performed.',flush=True)


if __name__=='__main__':main()
