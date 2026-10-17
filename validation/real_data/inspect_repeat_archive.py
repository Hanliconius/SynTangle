"""Inventory deposited repeat ZIP members without silently choosing a track."""
import argparse,json,zipfile
from pathlib import Path
from prepare_published_inputs import download
from prepare_rhynchospora_inputs import verify

def main():
    p=argparse.ArgumentParser();p.add_argument('--prior',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    records=json.loads((a.prior/'repeat_candidates.json').read_text())
    records=[r for r in records if r['id']==333740]
    if len(records)!=1 or records[0]['restricted']:raise ValueError('Expected public deposited repeat ZIP')
    record=records[0];a.output.mkdir(parents=True,exist_ok=True);archive=a.output/'repeat_annotation.zip'
    download('https://edmond.mpg.de/api/access/datafile/333740',archive);verify(archive,record)
    members=[]
    with zipfile.ZipFile(archive) as z:
        # CRC validation streams members; no bulk extraction.
        bad=z.testzip()
        if bad:raise ValueError('Corrupt ZIP member: '+bad)
        for info in z.infolist():
            if info.is_dir():continue
            row=dict(member=info.filename,bytes=info.file_size)
            if info.filename.lower().endswith(('.gff','.gff3','.bed','.out','.txt')):
                with z.open(info) as f:
                    row['preview']=[f.readline().decode('utf-8',errors='replace').rstrip()[:400] for _ in range(4)]
            members.append(row)
    (a.output/'archive_inventory.json').write_text(json.dumps(dict(source=record,members=members),indent=2)+'\n')
    candidates=[r for r in members if 'tenuis' in r['member'].lower() and r['member'].lower().endswith(('.gff','.gff3','.bed','.gff.gz','.gff3.gz','.bed.gz'))]
    for r in candidates[:12]:print('TRACK',r['member'],'bytes='+str(r['bytes']),flush=True)
    print('INVENTORY',a.output/'archive_inventory.json','members='+str(len(members)),flush=True)
    (a.output/'COMPLETE').write_text('PASS archive inspection only; no filtering performed\n')
if __name__=='__main__':main()
