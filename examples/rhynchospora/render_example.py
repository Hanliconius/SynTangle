"""Render the archived, independently rescored two comparison layouts; no optimization."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch,Rectangle
from syntangle.table_input import load_tables
from syntangle.saved_layout import decode_saved_layout,validate_saved_layout
from syntangle.orientation_space import orientation_basis
from syntangle.public_plot import directed_endpoints
from syntangle.visualize import _unambiguous_links
from syntangle.workflow import metrics

base=Path(__file__).resolve().parent
fixture,_=load_tables(base/'chromosomes.tsv',base/'links.tsv')
audit=json.loads((base/'reference_layouts.json').read_text())
cs={c.ref:c for c in fixture.chromosomes}
basis=orientation_basis(fixture,fixture.chromosome_refs)
p=audit['provenance'];palette=p['reference_palette'];colorref=p['homology_colour_reference']
maxlen=max(sum(c.length for c in fixture.chromosomes if c.ref.species_id==sp) for sp in fixture.species_ids)
gap=.015*maxlen
width=max(sum(c.length+gap for c in fixture.chromosomes if c.ref.species_id==sp) for sp in fixture.species_ids)
links=[(i,link) for i,(a,b) in enumerate(zip(fixture.species_ids,fixture.species_ids[1:])) for link in _unambiguous_links(fixture,a,b)]
assert len(links)==37
plt.rcParams['svg.fonttype']='none'
fig,axes=plt.subplots(1,2,figsize=(14,4.8))
titles=['Initial sequence order','SynTangle optimized']
panels=[audit['panels'][0],audit['panels'][3]]
for ax,panel,title in zip(axes.flat,panels,titles):
 state=decode_saved_layout(fixture,panel['state']);validate_saved_layout(fixture,state,basis)
 score=metrics(fixture,state);assert score['crossings']==panel['crossings']
 offsets={}
 for i,sp in enumerate(fixture.species_ids):
  x=0
  for ref in state.chromosome_order[sp]:
   c=cs[ref];offsets[ref]=x
   ax.add_patch(Rectangle((x,-i-.027),c.length,.054,color='#27333D',zorder=3))
   label=ref.chromosome_id.replace('scaffold_','').replace('Chr','')
   ax.text(x+c.length/2,-i+.09,label+(' +' if state.chromosome_orientation[ref]>0 else ' −'),ha='center',fontsize=8)
   x+=c.length+gap
 for i,(h,(c1,b1),(c2,b2)) in links:
  a,b=directed_endpoints(c1,b1,state.chromosome_orientation[c1.ref]);c,d=directed_endpoints(c2,b2,state.chromosome_orientation[c2.ref])
  a+=offsets[c1.ref];b+=offsets[c1.ref];c+=offsets[c2.ref];d+=offsets[c2.ref]
  y,z=-i,-i-1;mid=(y+z)/2
  vertices=[(a,y),(a,mid),(c,mid),(c,z),(d,z),(d,mid),(b,mid),(b,y),(a,y)]
  codes=[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY]
  ax.add_patch(PathPatch(MPath(vertices,codes),facecolor=palette[colorref[h]],alpha=.65,edgecolor='none'))
 ax.set_title(f'{title}\n{score["crossings"]} midpoint crossings',loc='left',fontsize=12)
 ax.set_xlim(-width*.02,width);ax.set_ylim(-2.17,.37)
 ax.set_yticks([0,-1,-2],['R. tenuis','R. austrobrasiliensis','R. breviuscula'],fontsize=9,fontstyle='italic')
 ax.set_xticks([]);ax.tick_params(length=0)
 for spine in ax.spines.values():spine.set_visible(False)
fig.suptitle('Rhynchospora: fewer crossings, identical evidence',fontsize=16)
fig.text(.5,.025,'55 → 9 crossings (83.6% reduction) · numerical lower bound 9 · crossing gap 0\n37/37 links and 10/10 chromosomes retained · whole-chromosome moves only\nCrossing optimality applies to midpoint crossings, not ribbon twisting.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.17,1,.93))
out=base.parents[1]/'docs/images';out.mkdir(parents=True,exist_ok=True)
for ext in ['svg','pdf']:fig.savefig(out/f'rhynchospora-example.{ext}')
plt.close(fig)
print('Verified and rendered archived baseline and optimized layouts: 55, 9 crossings')
