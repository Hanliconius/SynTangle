"""Recover the paper's chromosome-sized scaffolds and reuse completed Liftoff."""
import argparse,json,re
from pathlib import Path
from types import SimpleNamespace
from leptidea_liftoff import LEPTIDEA,sha,save,collect
# Male assemblies, Table 1, DOI 10.1007/s10577-023-09713-z.
TABLE={'LsinapisSpaM':(52,16290,.904),'LsinapisSweM':(29,16022,.902),
       'LjuvernicaM':(42,16073,.907),'LrealiM':(26,15448,.884)}

def select(rows,name):
    count,total,fraction=TABLE[name]
    if len(rows)!=total:raise ValueError('Scaffold count differs from paper Table 1: '+name)
    if any(r.get('assemblyAccession')!=LEPTIDEA[name] for r in rows):raise ValueError('Wrong assembly')
    numbered={}
    for r in rows:
        m=re.fullmatch(r'HiC_scaffold_(\d+)',r.get('sequenceName',''))
        if m:
            number=int(m[1])
            if number in numbered:raise ValueError('Duplicate scaffold name')
            numbered[number]=r
    if not set(range(1,count+1))<=set(numbered):raise ValueError('Paper-sized numbered scaffold series absent: '+name)
    chosen=[numbered[i] for i in range(1,count+1)]
    lengths=[int(r['length']) for r in chosen]
    actual=sum(lengths)/sum(int(r['length']) for r in rows)
    if abs(actual-fraction)>.0015:raise ValueError(f'Pseudochromosome fraction differs: {name} {actual:.5f} vs {fraction}')
    if any(int(r['length'])>min(lengths) for r in rows if r not in chosen):raise ValueError('Unselected scaffold exceeds smallest selected pseudochromosome')
    chromosomes=[dict(accession=r['genbankAccession'],sequence_name=r['sequenceName'],length=int(r['length'])) for r in chosen]
    if len({r['accession'] for r in chromosomes})!=count:raise ValueError('Duplicate accession')
    return dict(accession=LEPTIDEA[name],chromosomes=chromosomes,total_scaffolds=len(rows),selected_fraction=actual,published_fraction=fraction,numbered_lengths_descending=lengths==sorted(lengths,reverse=True))

def main():
    p=argparse.ArgumentParser();p.add_argument('--prior',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    mapping={'source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC9859909/#Tab1',
      'author_sort_script':'https://github.com/EBC-butterfly-genomics-team/Leptidea_chromosome_research2022/blob/8612edcad6d79d70c6eba591f52be36db1b89a02/Scripts/Assembly/finalize_assembly.sh',
      'scope':'Recovered chromosome-sized scaffold set; HiC labels retained. Paper figure order and original synteny ribbons are not established.', 'assemblies':{}}
    for name in LEPTIDEA:
        report=a.prior/'targets'/name/'sequence_report.jsonl'
        selected=select([json.loads(l) for l in report.read_text().splitlines() if l.strip()],name)
        selected['sequence_report_sha256']=sha(report);mapping['assemblies'][name]=selected
        print(f"CHROMOSOMES {name}: {len(selected['chromosomes'])}; fraction={selected['selected_fraction']:.4f}",flush=True)
    save(a.output/'chromosome_map.json',mapping)
    for name in ('targets','annotations','donor'):
        (a.output/name).symlink_to((a.prior/name).resolve(),target_is_directory=True)
    collect(SimpleNamespace(run=a.output,chromosome_map=a.output/'chromosome_map.json'))
    (a.output/'COMPLETE').write_text('PASS chromosome recovery and candidate homology collection\n')
if __name__=='__main__':main()
