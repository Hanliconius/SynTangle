"""Render actual saved layouts; do not optimize or infer missing benchmark states."""
import argparse,csv,json,re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,PathPatch
from matplotlib.path import Path as MPath
from matplotlib.backends.backend_pdf import PdfPages
from syntangle import load_validation_bundle,score_crossings
from syntangle.layout import LayoutState,initial_layout_state
from syntangle.saved_layout import decode_saved_layout,validate_saved_layout
from syntangle.orientation_space import orientation_basis
from syntangle.visualize import _unambiguous_links
from syntangle.model import Fixture,Chromosome,ChromosomeRef,BlockOccurrence

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
COLORS=list(plt.get_cmap('tab20').colors)

def colour(h):
    m=re.search(r'(?:gene|H)(\d+)',h)
    return COLORS[(int(m.group(1))-1 if m else sum(map(ord,h)))%20]


def chromosome_geometry(fixture,state,*,independent_rows=False):
    """Display coordinates only; leave evidence and layout states untouched."""
    chroms={c.ref:c for c in fixture.chromosomes};gap=max(c.length for c in fixture.chromosomes)*.08
    spans={sp:sum(chroms[r].length for r in state.chromosome_order[sp])+gap*(len(state.chromosome_order[sp])-1) for sp in fixture.species_ids}
    common_span=max(spans.values())
    geom={}
    for i,sp in enumerate(fixture.species_ids):
        span=spans[sp] if independent_rows else common_span
        x=0
        for ref in state.chromosome_order[sp]:
            geom[ref]=(x/span,chroms[ref].length/span,i);x+=chroms[ref].length+gap
    return geom


def riparian(ax,fixture,state,title,subtitle,*,independent_rows=False,respect_block_strand=False):
    validate_saved_layout(fixture,state,orientation_basis(fixture,fixture.chromosome_refs))
    score=score_crossings(fixture,state).crossings
    chroms={c.ref:c for c in fixture.chromosomes}
    geom=chromosome_geometry(fixture,state,independent_rows=independent_rows)
    def interval(c,b):
        x,w,y=geom[c.ref];a=b.start/c.length;d=b.end/c.length
        if state.chromosome_orientation[c.ref]<0:a,d=1-d,1-a
        return x+a*w,x+d*w,y
    for a,b in zip(fixture.species_ids,fixture.species_ids[1:]):
        for h,(c1,b1),(c2,b2) in _unambiguous_links(fixture,a,b):
            l1,r1,y1=interval(c1,b1);l2,r2,y2=interval(c2,b2)
            if respect_block_strand and b1.strand == "-":l1,r1=r1,l1
            if respect_block_strand and b2.strand == "-":l2,r2=r2,l2
            mid=(y1+y2)/2
            vertices=[(l1,y1+.065),(l1,mid),(l2,mid),(l2,y2-.065),(r2,y2-.065),(r2,mid),(r1,mid),(r1,y1+.065),(l1,y1+.065)]
            codes=[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY]
            ax.add_patch(PathPatch(MPath(vertices,codes),facecolor=colour(h),edgecolor=colour(h),linewidth=.12,alpha=.23))
    for ref,(x,w,y) in geom.items():
        ax.add_patch(Rectangle((x,y-.065),w,.13,facecolor='#eff1f3',edgecolor='#384553',linewidth=.35,zorder=3))
        for b in chroms[ref].blocks:
            l,r,_=interval(chroms[ref],b)
            ax.add_patch(Rectangle((l,y-.053),r-l,.106,facecolor=colour(b.homology_id),edgecolor='none',zorder=4))
    ax.set_xlim(-.055,1.005);ax.set_ylim(len(fixture.species_ids)-.65,-.55)
    ax.set_yticks(range(len(fixture.species_ids)),[sp.replace('species','') for sp in fixture.species_ids]);ax.tick_params(axis='y',length=0,labelsize=8)
    ax.set_xticks([])
    for s in ax.spines.values():s.set_visible(False)
    ax.set_title(title+'\n'+f'C = {score:,}  |  '+subtitle,loc='left',fontsize=10,pad=10)
    return score


def verify(fixture,data,expected=None):
    state=decode_saved_layout(fixture,data)
    validate_saved_layout(fixture,state,orientation_basis(fixture,fixture.chromosome_refs))
    count=score_crossings(fixture,state).crossings
    if expected is not None and count!=int(expected):raise AssertionError((count,expected))
    return state


def demonstration():
    chroms=[]
    for sp in ('speciesA','speciesB','speciesC'):
        for j in range(12):
            ref=ChromosomeRef(sp,str(j+1));blocks=tuple(BlockOccurrence(f'{sp}_{j}_{k}',f'gene{j+1}_{k}',15+k*20,23+k*20,'+') for k in range(3))
            chroms.append(Chromosome(ref,100,11-j if sp=='speciesB' else j,1,blocks))
    fixture=Fixture(1,'zero_demo','zero_demo','Synthetic paired presentation',tuple(chroms),(),{})
    after=LayoutState({sp:tuple(c.ref for c in chroms if c.ref.species_id==sp) for sp in fixture.species_ids},dict.fromkeys(fixture.chromosome_refs,1))
    return fixture,initial_layout_state(fixture),after


def main():
    p=argparse.ArgumentParser();p.add_argument('--data-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    inputs=args.data_root/'SynTangle_pegasus_test/local_results/genespace_stress'
    runs=list((args.data_root/'SynTangle_pipeline_test/local_results').glob('policy_*'))
    if len(runs)!=1:raise ValueError('Expected one policy run in supplied bundle')
    run=runs[0];audit=[];results=[]
    pdf=args.output/'SynTangle_visual_summary.pdf'
    with PdfPages(pdf) as report:
        f,before,after=demonstration();fig,axs=plt.subplots(1,2,figsize=(15,7))
        c0=riparian(axs[0],f,before,'A tangled presentation','same simulated biology')
        c1=riparian(axs[1],f,after,'A legal rearrangement of the drawing','C*=0, proved by nonnegativity')
        assert c1==0
        fig.suptitle('Visual tangledness can be entirely avoidable',fontsize=20,x=.055,ha='left')
        fig.text(.055,.075,f'{c0:,} crossings removed without changing any within-chromosome coordinates or homology links.',fontsize=12)
        fig.text(.055,.035,'Separate synthetic demonstration: 3 species, 12 chromosomes each. Colours group homologous links; bars are whole chromosomes.',fontsize=9,color='#546171')
        fig.subplots_adjust(left=.055,right=.975,top=.82,bottom=.16,wspace=.14)
        report.savefig(fig);fig.savefig(args.output/'zero_crossing_preview.png',dpi=120);plt.close(fig)
        for n,(case,level) in enumerate([('stress_6sp_20chr_5ev_random','Easy stress case'),('stress_8sp_30chr_8ev_random','Medium stress case'),('stress_10sp_40chr_12ev_random','Hard stress case')]):
            f=load_validation_bundle(inputs/case);comparison=inputs/'comparison'/case
            rows={r['method']:r for r in csv.DictReader((comparison/'results.tsv').open(),delimiter='\t')}
            panels=[]
            for method,label in [('input','Original input'),('GENESPACE','Native GENESPACE'),('GENESPACE_plus_flips','GENESPACE + our flip assistance')]:
                state=verify(f,json.loads((comparison/(method+'.layout.json')).read_text()),rows[method]['crossings'])
                panels.append((state,label,'display' if method=='input' else 'heuristic; multi-reference pilot'))
            short=json.loads((run/case/'mirror_selective_public_150_r0/result.json').read_text())
            panels.append((verify(f,short['optimized_state'],short['upper_bound']),'SynTangle - short allowance',f"{short['seconds']:.2f} s; "+('matching global bounds' if short['lower_bound']==short['upper_bound'] else f"L={short['lower_bound']:,}, U={short['upper_bound']:,}")))
            longfile=run/case/'mirror_selective_public_600_r0/result.json'
            if longfile.exists():
                long=json.loads(longfile.read_text());panels.append((verify(f,long['optimized_state'],long['upper_bound']),'SynTangle - long allowance',f"{long['seconds']:.2f} s; matching global bounds" if long['lower_bound']==long['upper_bound'] else 'unresolved'))
            else:long=short
            cols=2;nr=3 if n==2 else 2
            fig,axes=plt.subplots(nr,cols,figsize=(17,13 if nr==3 else 10));axes=list(axes.flat)
            counts=[]
            for ax,(state,label,note) in zip(axes,panels):counts.append(riparian(ax,f,state,label,note))
            for ax in axes[len(panels):]:
                ax.axis('off');ax.text(.05,.8,'The same evidence in every panel',fontsize=14,transform=ax.transAxes)
                ax.text(.05,.63,'Only whole-chromosome order and orientation change.\nAll homologous links are drawn; none are filtered.\nColours and chromosome scales are fixed within the case.',fontsize=10,transform=ax.transAxes,linespacing=1.6,va='top')
            title=f'{level}: {len(f.species_ids)} species | '+case.split('_')[2].replace('chr',' ancestral chromosomes')
            fig.suptitle(title,x=.055,ha='left',fontsize=20)
            fig.text(.055,.034,'Synthetic benchmark, random presentation. C counts adjacent-layer link-midpoint crossings. Palette repeats beyond 20 homology groups.\nAll rows may move in SynTangle; native GENESPACE uses reference-based variants. Difficulty labels describe these fixtures, not a general graph classification.',fontsize=9,color='#546171')
            fig.subplots_adjust(left=.055,right=.975,top=.9,bottom=.095,hspace=.32,wspace=.13)
            report.savefig(fig);fig.savefig(args.output/(('easy','medium','hard')[n]+'_riparian_preview.png'),dpi=110);plt.close(fig)
            results.append(dict(case=case,input=counts[0],genespace=counts[1],genespace_plus_flips=counts[2],short=short['upper_bound'],short_lower=short['lower_bound'],final=long['upper_bound'],final_lower=long['lower_bound'],solver_seconds=long['seconds']))
            for state,label,note in panels:audit.append(dict(case=case,panel=label,crossings=score_crossings(f,state).crossings,legal=True,link_count=score_crossings(f,state).link_count))
        fig,axs=plt.subplots(1,2,figsize=(15,8));methods=['input','genespace','genespace_plus_flips','short','final'];labels=['Input','GENESPACE','GS + flips','SynTangle short','SynTangle final']
        for i,row in enumerate(results):
            axs[0].plot(range(5),[row[k] for k in methods],marker='o',linewidth=2,label=row['case'].split('_')[1],color=['#0072B2','#D55E00','#009E73'][i])
        axs[0].set_yscale('log');axs[0].set_xticks(range(5),labels,rotation=25,ha='right');axs[0].set_ylabel('Crossings (log scale)');axs[0].legend(title='Species');axs[0].set_title('Crossing reduction on saved random presentations',loc='left',fontsize=12);axs[0].grid(axis='y',alpha=.2)
        overview=[]
        for line in (run/'summary.md').read_text().splitlines():
            if line.startswith('| 600 | public |'):
                v=[x.strip() for x in line.split('|')[1:-1]];overview.append((v[2],float(v[5])))
        if overview:
            names,values=zip(*overview);axs[1].barh(range(len(names)),values,color=['#0072B2' if 'mirror' in name else '#949daa' for name in names]);axs[1].set_yticks(range(len(names)),names);axs[1].invert_yaxis();axs[1].set_xlabel('Median worker wall time (seconds)');axs[1].set_title('Largest cases: public start, 600 s allowance',loc='left',fontsize=12)
        fig.suptitle('Layout quality and time to finish',x=.06,ha='left',fontsize=20)
        fig.text(.06,.07,'Left: canonically rescored saved layouts; short allowance 150 s, long allowance 600 s.\nRight: six runs per variant across three presentations, two repeats. A completed run need not have proved optimality.\nRuntime comparisons here are between SynTangle variants; historical GENESPACE timing scopes differ.',fontsize=10,color='#546171')
        fig.subplots_adjust(left=.07,right=.97,top=.86,bottom=.24,wspace=.6);report.savefig(fig);fig.savefig(args.output/'benchmark_summary_preview.png',dpi=120);plt.close(fig)
    (args.output/'layout_verification.json').write_text(json.dumps(dict(panels=audit,summary=results,zero_demo=dict(input=c0,final=c1),scope='canonical rescore and legal-state checks; historical audit verdicts are not rerun'),indent=2)+'\n')
    print(json.dumps(results,indent=2));print(pdf)

if __name__=='__main__':main()

