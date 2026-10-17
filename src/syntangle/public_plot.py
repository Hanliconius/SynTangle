"""Vector comparison drawing with immutable directed ribbon endpoints."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
from .visualize import _unambiguous_links


def directed_endpoints(chromosome, block, orientation):
    a,b = (block.start,block.end) if block.strand=='+' else (block.end,block.start)
    return (a,b) if orientation==1 else (chromosome.length-a,chromosome.length-b)


def plot_comparison(fixture, before, after, path):
    from .workflow import metrics
    chromosomes = {c.ref:c for c in fixture.chromosomes}
    max_length = max(sum(c.length for c in fixture.chromosomes if c.ref.species_id==sp) for sp in fixture.species_ids)
    gap = max_length*.015
    fig, axes = plt.subplots(1,2,figsize=(14,max(4,len(fixture.species_ids)*.75)),squeeze=False)
    links = [(i,link) for i,(a,b) in enumerate(zip(fixture.species_ids,fixture.species_ids[1:])) for link in _unambiguous_links(fixture,a,b)]
    colors = {h:plt.get_cmap('tab20')(i%20) for i,h in enumerate(dict.fromkeys(link[0] for _,link in links))}
    width = max(sum(chromosomes[r].length+gap for r in order) for order in before.chromosome_order.values())
    for ax,state,title in zip(axes[0],(before,after),('Input','Optimized')):
        offsets = {}
        for i,sp in enumerate(fixture.species_ids):
            x=0
            for ref in state.chromosome_order[sp]:
                offsets[ref]=x;c=chromosomes[ref]
                ax.add_patch(Rectangle((x,-i-.025),c.length,.05,facecolor='#444444',zorder=3))
                ax.text(x+c.length/2,-i+.08,ref.chromosome_id+(' +' if state.chromosome_orientation[ref]>0 else ' −'),ha='center',fontsize=7)
                x+=c.length+gap
        for i,(h,(c1,b1),(c2,b2)) in links:
            a,b=directed_endpoints(c1,b1,state.chromosome_orientation[c1.ref])
            c,d=directed_endpoints(c2,b2,state.chromosome_orientation[c2.ref])
            a+=offsets[c1.ref];b+=offsets[c1.ref];c+=offsets[c2.ref];d+=offsets[c2.ref]
            y,z=-i,-i-1;mid=(y+z)/2
            vertices=[(a,y),(a,mid),(c,mid),(c,z),(d,z),(d,mid),(b,mid),(b,y),(a,y)]
            codes=[Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.LINETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.CLOSEPOLY]
            ax.add_patch(PathPatch(Path(vertices,codes),facecolor=colors[h],edgecolor='none',alpha=.5))
        score=metrics(fixture,state)
        ax.set_title(f'{title}\n{score["crossings"]} midpoint crossings; {score["twists"]} twisted ribbons',fontsize=11)
        ax.set_xlim(-width*.02,width);ax.set_ylim(-len(fixture.species_ids)+.55,.4)
        ax.set_yticks([-i for i in range(len(fixture.species_ids))],fixture.species_ids)
        ax.tick_params(axis='both',length=0);ax.set_xticks([])
        for spine in ax.spines.values():spine.set_visible(False)
    fig.text(.5,.015,'Common coordinate scale · identical links and internal coordinates · whole-chromosome moves only',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.04,1,1));fig.savefig(path);plt.close(fig)
