"""Author-style Fig. 4f reconstruction from saved evidence; no solver runs."""
import argparse,hashlib,json,math,shutil,zipfile
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
from pypdf import PdfReader,PdfWriter,Transformation
from syntangle.fixtures import load_fixture
from syntangle.saved_layout import decode_saved_layout,validate_saved_layout
from syntangle.orientation_space import orientation_basis
from syntangle.layout import initial_layout_state,score_crossings
from syntangle.visualize import _unambiguous_links
from prepare_published_inputs import download

AUTHOR_COMMIT='91ae879d1d72386d393f38b5e5b4d111a2f1458b'
PDF_URL='https://www.mpipz.mpg.de/5631140/Ivankovic_Nat_Comm_2024.pdf'
SCRIPT_URL=f'https://raw.githubusercontent.com/Jeremias-Brand/PlanarianGenomeAnalysis/{AUTHOR_COMMIT}/scripts/fig4_busco.R'
NAMES=['schMedS3_h1','schMedS3_h2','schPol2','schNov1','schLug1','cloSin','schMan','hymMic','taeMul']
LABELS=['haplotype 1','haplotype 2','S. polychroa','S. nova','S. lugubris','C. sinensis','S. mansoni','H. microstoma','T. multiceps']
# R first maps chromosome -> hex string; scale_color_manual then remaps the
# alphabetically ordered hex-string factor levels onto its unnamed palette.
BASE=['#000000','#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7']
MANUAL=['#E69F00','#F0E442','#D55E00','#0072B2','#CC79A7','#000000','#009E73','#56B4E9']
REMAP=dict(zip(sorted(BASE),MANUAL))
COLOURS={c:REMAP[v] for c,v in zip([f'chr{i}' for i in range(1,8)]+['chrZW'],BASE)}

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def curve(x1,y1,x2,y2):
    # Exact control points in generate_bezier_curve(), rather than the previous
    # renderer's vertical-tangent ribbon approximation.
    curvature=min(abs(x2-x1)/abs(y2-y1),.3)
    mid=(y1+y2)/2;amplitude=(y2-y1)/2
    return [(x1,y1),(x1+(x2-x1)/10,mid+amplitude*curvature),
            (x1+9*(x2-x1)/10,mid-amplitude*curvature),(x2,y2)]

def geometry(fixture,state):
    chroms={c.ref:c for c in fixture.chromosomes};g={}
    for row,sp in enumerate(fixture.species_ids):
        for j,ref in enumerate(state.chromosome_order[sp]):g[ref]=(j*1.3+.3,row)
    def pos(c,b):
        x,y=g[c.ref];fraction=(b.start+b.end)/(2*c.length)
        if state.chromosome_orientation[c.ref]<0:fraction=1-fraction
        return x+fraction,y
    return g,chroms,pos

def draw(ax,fixture,state,provenance,labels=False):
    g,chroms,pos=geometry(fixture,state)
    colours=provenance['homology_colour_reference']
    # Grey chromosome bars are behind gene marks and coloured connecting lines.
    for ref,(x,y) in g.items():ax.plot([x,x+1],[y,y],color='#4D4D4D',linewidth=5*72/25.4,solid_capstyle='round',zorder=1)
    for a,b in zip(fixture.species_ids,fixture.species_ids[1:]):
        for h,(c,u),(d,v) in _unambiguous_links(fixture,a,b):
            x1,y1=pos(c,u);x2,y2=pos(d,v)
            path=MPath(curve(x1,y1,x2,y2),[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4])
            ax.add_patch(PathPatch(path,fill=False,edgecolor=COLOURS[colours[h]],linewidth=.5*72/25.4,alpha=.8,zorder=2))
    for ref,(x,y) in g.items():
        for b in chroms[ref].blocks:
            xx,_=pos(chroms[ref],b)
            ax.plot([xx,xx],[y-.07,y+.07],color=COLOURS[colours[b.homology_id]],linewidth=.5*72/25.4,zorder=3)
        if labels:ax.text(x+.5,y-.16,ref.chromosome_id+(' (-)' if state.chromosome_orientation[ref]<0 else ''),ha='center',fontsize=8)
    # Reconstruct the editorial phylogeny as topology only, with no branch lengths.
    ax.axvspan(12.05,13.25,ymin=.47,ymax=1,color='#E0E0E0',zorder=0)
    ax.axvspan(12.05,13.25,ymin=0,ymax=.41,color='#AAA6A6',zorder=0)
    def branch(left,right,x):
        for y in (left,right):ax.plot([x-.35,x],[y,y],color='#292929',lw=2)
        ax.plot([x,x],[left,right],color='#292929',lw=2)
        return (left+right)/2
    m=branch(0,1,12.25);m=branch(m,2,12.6)
    n=branch(3,4,12.6);branch(m,n,12.95)
    p=branch(5,6,12.25);q=branch(7,8,12.25);branch(p,q,12.6)
    # Group-joining backbone drawn separately; display tree carries no inference.
    ax.plot([12.95,13.4,13.4,12.6],[2.375,2.375,6.5,6.5],color='#292929',lw=2)
    for row,label in enumerate(LABELS):ax.text(11.9,row,label,ha='right',va='center',fontsize=10,fontstyle='italic' if row>1 else 'normal')
    ax.text(11.9,.5,'S. mediterranea',ha='right',va='center',fontsize=10,fontstyle='italic')
    ax.annotate('',xy=(5.28,4),xytext=(10,4),arrowprops=dict(arrowstyle='-|>',color='#D55E00',lw=2))
    ax.plot([10,10],[4,5],color='#D55E00',lw=2)
    ax.annotate('',xy=(9.18,5),xytext=(10,5),arrowprops=dict(arrowstyle='-|>',color='#D55E00',lw=2))
    ax.set_xlim(-.25,13.7);ax.set_ylim(8.3,-.35);ax.axis('off')

def write_panel(fixture,state,provenance,path,labels=False):
    fig,ax=plt.subplots(figsize=(10,10));draw(ax,fixture,state,provenance,labels)
    fig.subplots_adjust(left=.02,right=.98,top=.98,bottom=.02)
    fig.savefig(path,format='pdf');plt.close(fig)

def crop_published(source,out):
    page=PdfReader(source).pages[8]
    # Verified Fig. 4f bounds on journal PDF page 9, in PDF points.
    box=(245,225,548,502)
    page.cropbox.lower_left=box[:2];page.cropbox.upper_right=box[2:]
    page.mediabox.lower_left=box[:2];page.mediabox.upper_right=box[2:]
    w=PdfWriter();w.add_page(page);w.write(out)

def side_by_side(left,right,out):
    w=PdfWriter();canvas=w.add_blank_page(width=1100,height=720)
    for index,path in enumerate((left,right)):
        page=PdfReader(path).pages[0];box=page.mediabox
        x,y=float(box.left),float(box.bottom);width,height=float(box.width),float(box.height)
        scale=min(530/width,680/height)
        transform=Transformation().translate(-x,-y).scale(scale).translate(10+550*index,20)
        canvas.merge_transformed_page(page,transform,expand=False)
    w.write(out)

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    prepared=a.source/'prepared/flatworms_fig4f';fixture=load_fixture(prepared/'fixture.json')
    provenance=json.loads((prepared/'provenance.json').read_text())
    audit=json.loads((a.source/'figures/flatworms_fig4f/comparison_audit.json').read_text())
    if list(fixture.species_ids)!=NAMES or audit['fixture_sha256']!=digest(prepared/'fixture.json'):raise ValueError('Unexpected Fig. 4f evidence')
    basis=orientation_basis(fixture,fixture.chromosome_refs);original=initial_layout_state(fixture)
    for h,rows in provenance['original_coordinates'].items():
        for r in rows:
            c=next(c for c in fixture.chromosomes if c.ref.species_id==r['sp'] and c.ref.chromosome_id==r['chr'])
            b=next(b for b in c.blocks if b.homology_id==h)
            expected=round(r['start']/r['length'],3)
            if abs((b.start+b.end)/2-expected)>1.01e-6:raise ValueError('Coordinate differs from author rounding')
    panels=[]
    for i,item in enumerate(audit['panels']):
        state=decode_saved_layout(fixture,item['state']);validate_saved_layout(fixture,state,basis)
        score=score_crossings(fixture,state).crossings
        if score!=item['crossings']:raise ValueError('Saved crossings differ')
        if i==0 and state.to_dict()!=original.to_dict():raise ValueError('Saved first panel differs from native author order')
        output=a.output/f'panel_{i+1}.pdf';write_panel(fixture,state,provenance,output)
        panels.append({'title':item['title'],'crossings':score,'file':output.name,'state':item['state']})
        print('PANEL',i+1,'C',score,flush=True)
    write_panel(fixture,original,provenance,a.output/'reconstruction_labelled.pdf',True)
    # Also archive the exact publication crop, not a claim that rendering is identical.
    raw=a.output/'sources';raw.mkdir(exist_ok=True)
    download(PDF_URL,raw/'published.pdf');download(SCRIPT_URL,raw/'fig4_busco.R')
    crop_published(raw/'published.pdf',a.output/'published_panel4f.pdf')
    side_by_side(a.output/'published_panel4f.pdf',a.output/'panel_1.pdf',a.output/'published_vs_reconstructed.pdf')
    w=PdfWriter()
    for file in ('published_panel4f.pdf','panel_1.pdf','reconstruction_labelled.pdf','panel_4.pdf'):
        w.add_page(PdfReader(a.output/file).pages[0])
    w.write(a.output/'publication_and_reconstruction.pdf')
    report={'paper':'https://doi.org/10.1038/s41467-024-52380-9','author_script_commit':AUTHOR_COMMIT,
            'fixture_sha256':audit['fixture_sha256'],'published_pdf_sha256':digest(raw/'published.pdf'),
            'author_script_sha256':digest(raw/'fig4_busco.R'),'chromosome_gap':.3,'palette':COLOURS,
            'panels':panels,'retained_genes':len(provenance['original_coordinates']),
            'verified':['species order','native chromosome order and orientation','rounded gene-start coordinates','Schistosoma chromosome colours'],
            'limitations':['Visual similarity must be checked against archived publication crop; no pixel identity claim',
                          'Chromosome extents are equalised as in author code, not physical bp lengths',
                          'Phylogeny and red arrow are editorial reconstructions; branch lengths not inferred',
                          'Saved four-way comparison layouts are reused; no fresh solver execution']}
    (a.output/'reconstruction_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    (a.output/'ATTRIBUTION.txt').write_text('Fig. 4f reproduced from Ivankovic et al., Nature Communications 15, 8215 (2024), DOI 10.1038/s41467-024-52380-9. Original publication CC BY 4.0: https://creativecommons.org/licenses/by/4.0/. Reconstructed renderings and optimized layouts are adaptations; publication crop is the original. Author code: Jeremias-Brand/PlanarianGenomeAnalysis, commit '+AUTHOR_COMMIT+'.\n')
    archive=a.output.parent/'SynTangle_fig4f_reconstruction.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for f in sorted(a.output.rglob('*')):
            if f.is_file():z.write(f,str(f.relative_to(a.output)))
        z.write(Path(__file__),'code/'+Path(__file__).name)
    print('PDF:',a.output/'publication_and_reconstruction.pdf');print('ARCHIVE:',archive)

if __name__=='__main__':main()
