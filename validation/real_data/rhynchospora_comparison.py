"""Reconstruct transfer-derived Rhynchospora anchors on a cited pseudo-haplotype.

Not original published GENESPACE orthology or ribbon evidence. Gene ID equality
is meaningful only because both transfers use the exact same deposited donor.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from rhynchospora_followup import genes

SELECTED = ['scaffold_13_a1','scaffold_7_a1','scaffold_1_a1']
CASE = 'rhynchospora_transfer_reanalysis'
PAPER = 'https://www.nature.com/articles/s41586-026-11057-7'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def strict(path, donor_ids, lengths):
    records=list(genes(path)); counts=Counter(g['id'] for g in records)
    retained={}; excluded=Counter()
    for g in records:
        if g['sequence'] not in lengths or not 1<=g['start']<=g['end']<=lengths[g['sequence']]:
            raise ValueError('Transferred gene outside target sequence')
        attrs=g['attributes']
        if counts[g['id']]!=1 or g['id'] not in donor_ids:
            excluded['nonunique_or_unknown_donor']+=1;continue
        try: quality=[float(attrs[k]) for k in ('coverage','sequence_ID')]
        except (KeyError,ValueError):
            excluded['missing_quality']+=1;continue
        if not all(.9<=v<=1 for v in quality) or attrs.get('partial_mapping','False').lower()=='true':
            excluded['quality_or_partial']+=1;continue
        retained[g['id']]=g
    return retained,dict(excluded)


def prepare(a):
    from pyfaidx import Fasta
    prior=json.loads((a.prior/'manifest.json').read_text())
    inputs=Path(prior['inputs'])
    source=inputs/'sources/R_austrobrasiliensis_paper_deposit/Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta'
    if not (a.prior/'annotation/LIFTOFF_COMPLETE').exists():raise ValueError('Prior breviuscula transfer incomplete')
    out=a.run/'sources';out.mkdir()
    # Index a private symlink, never write an index alongside the deposited source.
    private=out/'all_austro.fasta';private.symlink_to(source)
    fasta=Fasta(str(private),indexname=str(out/'all_austro.fai'))
    if any(name not in fasta for name in SELECTED):raise ValueError('Cited scaffolds missing')
    with (out/'selected_austro.fasta').open('w') as f:
        for name in SELECTED:
            f.write('>'+name+'\n')
            for start in range(0,len(fasta[name]),1048560):
                sequence=str(fasta[name][start:start+1048560])
                f.write('\n'.join(sequence[i:i+80] for i in range(0,len(sequence),80))+'\n')
    fasta.close()
    provenance=dict(paper=PAPER,selection=SELECTED,
        selection_basis='Extended Data Fig. 4 caption explicitly uses scaffolds 1, 7, 13 as a pseudo-haplotype; NOT confirmation of exact Fig. 1a inputs',
        initial_order_basis='13,7,1 matches short-long-long bar geometry; orientations native and unverified',
        original_publication_evidence=False,prior_run=str(a.prior),inputs=str(inputs),
        selected_fasta_sha256=sha(out/'selected_austro.fasta'),
        donor_fasta_sha256=sha(a.prior/'sources/donor.fasta'),donor_gff_sha256=sha(a.prior/'sources/donor.gff'),
        existing_breviuscula_gff_sha256=sha(a.prior/'annotation/lifted.gff'))
    (a.run/'source_manifest.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('PREPARED cited austrobrasiliensis pseudo-haplotype: 13,7,1; previous breviuscula transfer reused')


def lift(a):
    out=a.run/'annotation';out.mkdir()
    cmd=['liftoff',str(a.run/'sources/selected_austro.fasta'),str(a.prior/'sources/donor.fasta'),
         '-g',str(a.prior/'sources/donor.gff'),'-o',str(out/'lifted.gff'),'-u',str(out/'unmapped.txt'),
         '-dir',str(out/'intermediate'),'-p',str(a.threads),'-a','0.5','-s','0.5','-exclude_partial']
    (out/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
    with (out/'liftoff.out').open('w') as stdout,(out/'liftoff.err').open('w') as stderr:
        subprocess.run(cmd,check=True,cwd=out,stdout=stdout,stderr=stderr)
    (out/'COMPLETE').write_text('PASS\n')
    print('COMPLETE austrobrasiliensis transfer; detailed logs retained in annotation/')


def fixture(a):
    from pyfaidx import Fasta
    from syntangle.fixtures import fixture_from_dict
    if not (a.run/'annotation/COMPLETE').exists():raise ValueError('Austrobrasiliensis transfer incomplete')
    # Private indices, including the previously validated source sequences.
    lengths={}
    for sp,path in [('R_tenuis',a.prior/'sources/donor.fasta'),
                    ('R_austrobrasiliensis',a.run/'sources/selected_austro.fasta'),
                    ('R_breviuscula',a.prior/'sources/target.fasta')]:
        f=Fasta(str(path),indexname=str(a.run/'sources'/f'{sp}.fai'))
        lengths[sp]={name:len(f[name]) for name in f.keys()};f.close()
    donor=list(genes(a.prior/'sources/donor.gff'))
    if len({g['id'] for g in donor})!=len(donor):raise ValueError('Duplicate donor IDs')
    donors={g['id']:g for g in donor}
    austro,ex_a=strict(a.run/'annotation/lifted.gff',donors,lengths['R_austrobrasiliensis'])
    brevi,ex_b=strict(a.prior/'annotation/lifted.gff',donors,lengths['R_breviuscula'])
    # All T-A anchors and all A-B shared-donor anchors are retained, with no sorting filter.
    scored=set(austro); triple=scored & set(brevi)
    if len(scored)<10 or len(triple)<10:raise ValueError('Too few strict transferred anchors for both adjacent pairs')
    mappings={'R_tenuis':donors,'R_austrobrasiliensis':austro,'R_breviuscula':brevi}
    data=dict(fixture_version=1,id=CASE,title='Rhynchospora reconstructed transfer-anchor comparison',
              purpose='Whole-chromosome ordering of explicit transfer-derived adjacent anchors',species=[],orientation_constraints=[])
    for sp,table in mappings.items():
        names=SELECTED if sp=='R_austrobrasiliensis' else list(lengths[sp])
        chromosomes=[]
        for rank,name in enumerate(names,1):
            blocks=[]
            for hid in sorted(scored):
                if hid not in table or table[hid]['sequence']!=name:continue
                g=table[hid]
                blocks.append(dict(occurrence_id=sp+'_'+hid,homology_id=hid,
                                   start=g['start']-1,end=g['end'],strand=g['strand']))
            chromosomes.append(dict(id=name,length=lengths[sp][name],display_rank=rank,display_orientation=1,blocks=blocks))
        data['species'].append(dict(id=sp,chromosomes=chromosomes))
    fixture_from_dict(data)
    provenance=json.loads((a.run/'source_manifest.json').read_text())
    provenance.update(source_kind='New Liftoff-derived gene anchors; not published ribbons',
        retained_links=len(scored)+len(triple),pair_audit=[dict(pair=['R_tenuis','R_austrobrasiliensis'],retained=len(scored)),
                                                       dict(pair=['R_austrobrasiliensis','R_breviuscula'],retained=len(triple))],
        orientation_note='Transfer quality audited; gene strands are not hard chromosome-orientation constraints.',
        chromosome_extent='Exact deposited FASTA lengths; all chromosomes in selected representations retained',
        colour_reference='R_tenuis',colour_note='Colours follow donor chromosome membership, constant across panels.',
        input_panel_title='Candidate paper order: reconstructed anchors',
        input_layout_note='Austro scaffolds 13,7,1: length-based candidate order; native orientations. NOT exact published drawing.',
        crossing_unit='gene-anchor',anchor_description='transfer-derived gene anchors',
        workflow_note='GENESPACE ordering functions on reconstructed anchors; NOT the original orthology/collinearity workflow.',
        objective='Unweighted strict gene-midpoint crossing count, fixed Tenuis-Austro-Brevi row chain',
        scope='Full selected chromosome representations; strict unique >=90% identity/coverage transfer anchors; no weak-hit filtering in solver',
        homology_warning='Shared donor IDs establish transfer correspondence, NOT independent reciprocal orthology validation.',
        strict_austro_genes=len(austro),strict_brevi_genes=len(brevi),triple_genes=len(triple),
        brevi_only_nonadjacent_genes=len(set(brevi)-set(austro)),transfer_exclusions={'austro':ex_a,'brevi':ex_b})
    out=a.run/'prepared'/CASE;out.mkdir(parents=True)
    (out/'fixture.json').write_text(json.dumps(data,indent=2)+'\n')
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    shutil.copyfile(out/'provenance.json',a.run/'anchor_audit.json')
    print(f'PREPARED {CASE}: T-A={len(scored)} A-B={len(triple)}; full selected chromosomes; transfer-derived correspondence only')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=['prepare','lift','fixture'],required=True)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--prior',type=Path,required=True)
    p.add_argument('--threads',type=int,default=8)
    a=p.parse_args();globals()[a.stage](a)
