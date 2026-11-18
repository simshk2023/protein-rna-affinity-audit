"""Wide (178 mm) recomposition of manuscript Figure 1.

Panels a and b are re-laid out from the frozen
figure1_bright_semicircle_20261101_v3 source: the native point sizes
(7.5/8/8.5/9/11 pt) and every text string are untouched; only axes rectangles
and text anchor coordinates change. Panel c is re-rendered at a larger canvas
by render_response_wide.py (same fonts/sizes/data) and composited with
PyMuPDF at 1:1 scale, cropped to its own ink bounding box so that the lower
right of the 178 mm page is used well.

Run:  python render_figure1_wide.py   [--reuse-response]
"""
from pathlib import Path
import json,sys,subprocess,hashlib,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap,Normalize
from scipy.stats import gaussian_kde
import pymupdf as fitz
from PIL import Image,ImageOps

O=Path(__file__).resolve().parent;R=O.parents[1]
DATA_SHA256='6dfeefd6cc8a1469fc227103ad2713b7135e7cd8d18bf02f5a5eba66dd866e8a'
if '--reuse-response' not in sys.argv:
    subprocess.run([sys.executable,str(O/'render_response_wide.py'),'--export'],check=True)
assert (O/'response_wide.pdf').exists(),'Run render_response_wide.py first.'
D=json.loads((O/'family_aggregates.json').read_text())
assert hashlib.sha256((O/'family_aggregates.json').read_bytes()).hexdigest()==DATA_SHA256,'Frozen family_aggregates.json changed.'
W,H=178,164
MM=72/25.4

# ---------------------------------------------------------------- panel c prep
# Measure the re-rendered panel-c ink extent, crop to it and decide the target
# rectangle before the base figure is built (the c heading is anchored to it).
resp=fitz.open(O/'response_wide.pdf');rp=resp[0]
rp_w,rp_h=rp.rect.width/MM,rp.rect.height/MM
assert abs(rp_w-104)<.01 and abs(rp_h-88)<.01,(rp_w,rp_h)
RDPI=600
pix=rp.get_pixmap(matrix=fitz.Matrix(RDPI/72,RDPI/72),alpha=False)
arr=np.frombuffer(pix.samples,dtype=np.uint8).reshape(pix.height,pix.width,pix.n)[:,:,:3]
ys,xs=np.where(arr.mean(axis=2)<250)
ink_bbox=(xs.min()/RDPI*25.4,ys.min()/RDPI*25.4,(xs.max()+1)/RDPI*25.4,(ys.max()+1)/RDPI*25.4)
resp_text=[s for b in rp.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
tx0=min(s['bbox'][0] for s in resp_text)/MM;ty0=min(s['bbox'][1] for s in resp_text)/MM
tx1=max(s['bbox'][2] for s in resp_text)/MM;ty1=max(s['bbox'][3] for s in resp_text)/MM
U=(min(ink_bbox[0],tx0),min(ink_bbox[1],ty0),max(ink_bbox[2],tx1),max(ink_bbox[3],ty1))
PAD=1.2
clip_mm=(U[0]-PAD,U[1]-PAD,U[2]+PAD,U[3]+PAD)
clip=fitz.Rect(*(v*MM for v in clip_mm))
c_w,c_h=clip_mm[2]-clip_mm[0],clip_mm[3]-clip_mm[1]
# Lower-right placement: ink right edge at C_RIGHT, ink top edge at C_TOP.
C_RIGHT,C_TOP=174.0,90.8
c_tx0=C_RIGHT-(U[2]-clip_mm[0])
c_ty0=C_TOP-(U[1]-clip_mm[1])
assert c_tx0>90 and c_ty0>80,(c_tx0,c_ty0)
assert c_tx0+c_w<=W-1 and c_ty0+c_h<=H-1,(c_tx0+c_w,c_ty0+c_h)
for s in resp_text:
    bb=fitz.Rect(s['bbox'])
    assert clip.contains(bb),'Panel-c crop would cut text: %r'%s['text']
target=fitz.Rect(c_tx0*MM,c_ty0*MM,(c_tx0+c_w)*MM,(c_ty0+c_h)*MM)
C_HEAD_X=c_tx0-1.2

# ---------------------------------------------------------------- base figure
from bright_palette import BLUE, TEAL, CORAL, AMBER, INK, SECONDARY, AXIS, GRID, LIGHT_BLUE, LIGHT_TEAL
GRAY=SECONDARY
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.linewidth':.55,'axes.edgecolor':AXIS,'xtick.major.size':2.5,'ytick.major.size':2.5,'xtick.major.width':.5,'ytick.major.width':.5,'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'stixsans','text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
def txt(x,y,s,fs=8,**kw):return fig.text(x/W,1-y/H,s,fontsize=fs,**kw)
def rect(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def heading(letter,x,y,s):
    txt(x,y,letter,11,weight='bold');txt(x+6.4,y,s,9,weight='bold')
def clean(ax):
    ax.spines[['top','right']].set_visible(False);ax.tick_params(pad=3)

# a. Three prediction-error columns widened to span the 178 mm canvas.
#    Native geometry (y=20..49.5 scatter, y=58..69 reduction) is preserved.
ROW2=74.5
COLS=[15.5,72.0,128.5];COLW=47.5
heading('a',4,4.5,'Prediction error')
fig.legend(handles=[Line2D([],[],ls='',marker='o',mfc='white',mec=TEAL,ms=4,label='Additive'),Line2D([],[],ls='',marker='^',mfc=BLUE,mec='white',mew=.3,ms=4.8,label='CoPRA'),Line2D([],[],ls='',marker='D',mfc=AMBER,mec='white',mew=.3,ms=4.5,label='Mean')],loc='upper right',bbox_to_anchor=(.99,.995),ncol=3,fontsize=7.5,frameon=False,handletextpad=.4,columnspacing=.9)
for i,(x,v,name) in enumerate(zip(COLS,D['B'],['Absolute MAE',r'$\Delta$MAE','cMAE'])):
    width=COLW;txt(x+width/2,12,name,8.5,ha='center',weight='bold')
    txt(x+width/2,16,f"{v['n']} {'units' if i==0 else 'pairs'} · {v['nf']} families",7.5,ha='center',color=GRAY)
    ax=rect(x,20,width,29.5);clean(ax)
    jitter=np.random.default_rng(33+i).uniform(-.16,.16,v['nf'])
    a=np.array(v['values']['baseline']);c=np.array(v['values']['copra'])
    if i==2:
        for k in range(len(a)):ax.plot([jitter[k],1+jitter[k]],[a[k],c[k]],color=LIGHT_BLUE,lw=.5,alpha=.65,zorder=0)
    ax.scatter(jitter,a,s=11,marker='o',facecolor='white',edgecolor=TEAL,lw=.6,alpha=.8)
    ax.scatter(1+jitter,c,s=15,marker='^',c=BLUE,edgecolor='white',lw=.25,alpha=.8)
    ax.scatter([0,1],[a.mean(),c.mean()],s=26,marker='D',c=AMBER,edgecolor='white',lw=.5,zorder=5)
    ax.set_xlim(-.4,1.4);ax.set_ylim(0,[5,3,1][i]);ax.set_xticks([0,1],['Additive','CoPRA']);ax.set_yticks([[0,1,2,3,4,5],[0,1,2,3],[0,.5,1]][i])
    if i==0:ax.set_ylabel('Family error (pKd)',labelpad=4)
    ae=rect(x,58,width,11);clean(ae);ae.axhline(0,color=GRAY,lw=.65,ls=(0,(3,2)))
    boot=np.array(v['bootstrap']);yy=np.linspace(boot.min(),boot.max(),240);den=gaussian_kde(boot)(yy);den=den/den.max()*.28
    ae.fill_betweenx(yy,.45,.45+den,color=BLUE,alpha=.22,lw=0)
    lo,hi=v['ci'];ae.plot([.43,.43],[lo,hi],color=BLUE,lw=1.5,zorder=3);ae.scatter(.43,v['gain'],marker='D',s=27,c=AMBER,edgecolor='white',lw=.5,zorder=4)
    lower=min(0,float(boot.min()));upper=max(0,float(boot.max()));pad=.08*(upper-lower)
    ae.set_xlim(0,1);ae.set_ylim(lower-pad,upper+pad);ae.set_xticks([]);ae.set_yticks([[0,.2],[0,.1],[-.05,0]][i]);ae.spines['bottom'].set_visible(False)
    if i==0:ae.set_ylabel('Reduction\n(pKd)',labelpad=4)

# b. Error-reduction matrix; wider 21.5 mm columns so the panel fills the
#    left part of the bottom row. Row height stays at the native 2.8 mm.
MAT_X,MAT_W=11.0,64.5          # 3 x 21.5 mm columns
MAT_TOP=ROW2+14.0              # 88.5
N_POS_X=MAT_X+MAT_W+6.5        # 82.0
ROW_H=2.8
heading('b',4,ROW2,'Error reduction');txt(10.4,ROW2+4.5,'220 pairs · 20 families',7.5,color=GRAY)
rows=sorted(D['local'],key=lambda r:-r['gain'][0]);G=np.array([r['gain'] for r in rows]);CMAX=abs(G).max(axis=0)
cmap=LinearSegmentedColormap.from_list('gain',[CORAL,'#FFFFFF',BLUE])
ax=rect(MAT_X,MAT_TOP,MAT_W,56);ax.set_xlim(0,3);ax.set_ylim(20,0);ax.set_xticks([]);ax.set_yticks(np.arange(20)+.5,[r['id'] for r in rows]);ax.tick_params(axis='y',length=0,pad=2.5,labelsize=7.5)
for s in ax.spines.values():s.set_visible(False)
def cell(ax,j,y,v):
    frac=v/CMAX[j];col=cmap(.5+.5*np.clip(frac,-1,1))
    ax.add_patch(Rectangle((j,y),1,1,fc=col,ec='white',lw=.8))
    ax.text(j+.5,y+.5,'+' if v>0 else '−',ha='center',va='center',fontsize=8,color=INK)
for i,row in enumerate(G):
    for j,v in enumerate(row):cell(ax,j,i,v)
for j,name in enumerate(['MAE',r'$\Delta$MAE','cMAE']):
    xc=MAT_X+(j+.5)*MAT_W/3;txt(xc,MAT_TOP-5.5,name,7.5,ha='center');txt(xc,MAT_TOP-2.5,f'±{CMAX[j]:.2f}',7.5,ha='center',color=GRAY)
txt(N_POS_X,MAT_TOP-5.5,r'$n_{pos}$',7.5,ha='center')
for i,r in enumerate(rows):txt(N_POS_X,MAT_TOP+(i+.5)*ROW_H,str(r['positions']),7.5,ha='center',va='center',color=GRAY)
MEAN_Y=MAT_TOP+58.2
am=rect(MAT_X,MEAN_Y,MAT_W,ROW_H);am.set_xlim(0,3);am.set_ylim(1,0);am.set_xticks([]);am.set_yticks([.5],['Mean']);am.tick_params(axis='y',length=0,pad=2.5,labelsize=7.5);am.get_yticklabels()[0].set_fontweight('bold')
for s in am.spines.values():s.set_visible(False)
for j,v in enumerate(G.mean(axis=0)):cell(am,j,0,v)
KEY_Y=MEAN_Y+5.0
key=rect(MAT_X,KEY_Y,MAT_W,1.4);key.pcolormesh(np.linspace(-1,1,257),[0,1],np.linspace(-1,1,256)[None,:],cmap=cmap,norm=Normalize(-1,1),antialiased=False,shading='flat',edgecolors='face',linewidth=.05);key.set_axis_off()
for x,s,ha in [(MAT_X,'−','left'),(MAT_X+MAT_W/2,'0','center'),(MAT_X+MAT_W,'+','right')]:txt(x,KEY_Y+4.2,s,7.5,ha=ha,color=GRAY)
txt(MAT_X+MAT_W/2,ROW2+85,'Column-scaled',7.5,ha='center',color=GRAY);txt(MAT_X+MAT_W/2,ROW2+88,'reduction (pKd)',7.5,ha='center',color=GRAY)

# c. Heading only; the plot itself is composited from the wide response PDF.
heading('c',C_HEAD_X,ROW2,'Amplitude and ordering')

sys.path.insert(0,str(Path(__file__).resolve().parent))
from visual_qa import audit_layout
issues=audit_layout(fig);assert not issues,issues
(O/'base_layout_wide_audit.json').write_text(json.dumps(issues,indent=2))
fig.savefig(O/'base_layout_wide.pdf');plt.close(fig)

# ---------------------------------------------------------------- composite
base=fitz.open(O/'base_layout_wide.pdf');page=base[0]
assert abs(page.rect.width-W*MM)<.01 and abs(page.rect.height-H*MM)<.01
page.show_pdf_page(target,resp,0,clip=clip)
base.save(O/'figure1.pdf',garbage=4,deflate=True);base.close();resp.close()

# matplotlib >=3.10 emits the non-BMP math-italic glyphs of `$n_{pos}$` as
# vector paths (visually identical) without a Unicode text layer.  Add an
# invisible render-mode-3 text copy with the exact Unicode code points and the
# native 7.5/5.25 pt sizes so that text extraction/accessibility matches the
# manuscript baseline.  The rendered pixels are asserted unchanged.
def add_npos_text_layer(pdf_path,anchor_x,anchor_y):
    import matplotlib as mpl
    fontfile=Path(mpl.get_data_path())/'fonts/ttf/STIXGeneralItalic.ttf'
    assert fontfile.exists(),fontfile
    doc=fitz.open(pdf_path);pg=doc[0]
    pix0=pg.get_pixmap(matrix=fitz.Matrix(200/72,200/72),alpha=False)
    fonts_before={f[3] for f in pg.get_fonts(full=True)}
    tw=fitz.TextWriter(pg.rect)
    font=fitz.Font(fontfile=str(fontfile))
    # Anchors follow the mathtext layout: main glyph then the 5.25 pt subscript.
    tw.append(fitz.Point((anchor_x-2.12)*MM,anchor_y*MM),'\U0001D62F',font=font,fontsize=7.5)
    tw.append(fitz.Point((anchor_x-0.83)*MM,(anchor_y+0.55)*MM),'pos',font=font,fontsize=5.25)
    tw.write_text(pg,render_mode=3)
    added=[f for f in pg.get_fonts(full=True) if f[3] not in fonts_before]
    assert len(added)==1,added
    m=re.search(r'/ToUnicode (\d+) 0 R',doc.xref_object(added[0][0]));assert m,'no ToUnicode stream'
    tufont=added[0][3].replace(' ','-')
    traces=[(t['size'],[c[1] for c in t['chars']]) for t in pg.get_texttrace() if t['font']==tufont]
    main=[t for t in traces if abs(t[0]-7.5)<.05][0];sub=[t for t in traces if abs(t[0]-5.25)<.05][0]
    assert len(main[1])==1 and len(sub[1])==3,(main,sub)
    target={main[1][0]:'\U0001D62F',sub[1][0]:'\U0001D631',sub[1][1]:'\U0001D630',sub[1][2]:'\U0001D634'}
    bfchar='\n'.join('<%04x> <%s>'%(g,ch.encode('utf-16-be').hex()) for g,ch in target.items())
    cmap=(b'/CIDInit /ProcSet findresource begin\n12 dict begin\nbegincmap\n'
          b'/CIDSystemInfo << /Registry (Adobe) /Ordering (UCS) /Supplement 0 >> def\n'
          b'/CMapName /Adobe-Identity-UCS def\n/CMapType 2 def\n1 begincodespacerange\n<0000> <ffff>\nendcodespacerange\n'
          +('%d beginbfchar\n%s\nendbfchar\n'%(len(target),bfchar)).encode('ascii')
          +b'endcmap\nCMapName currentdict /CMap defineresource pop\nend\nend\n')
    doc.update_stream(int(m.group(1)),cmap)
    updated = doc.tobytes(garbage=4, deflate=True)
    doc.close()
    pdf_path.write_bytes(updated)
    doc=fitz.open(pdf_path)
    pix1=doc[0].get_pixmap(matrix=fitz.Matrix(200/72,200/72),alpha=False)
    same=np.array_equal(np.frombuffer(pix0.samples,dtype=np.uint8),np.frombuffer(pix1.samples,dtype=np.uint8))
    doc.close()
    assert same,'invisible text layer changed the rendering'
    return ''.join(sorted(target.values()))

# Pre-layer vector SVG (visually identical to the final PDF).
pdf=fitz.open(O/'figure1.pdf')
(O/'figure1.svg').write_text(pdf[0].get_svg_image(text_as_path=False),encoding='utf-8')
pdf.close()
added_chars=add_npos_text_layer(O/'figure1.pdf',N_POS_X,MAT_TOP-5.5)

pdf=fitz.open(O/'figure1.pdf');page=pdf[0]
pw,ph=page.rect.width/MM,page.rect.height/MM
assert abs(pw-W)<.01 and abs(ph-H)<.01,(pw,ph)
assert not page.get_images(),'raster image in output'
fonts={r[3]:len(pdf.extract_font(r[0])[3]) for r in page.get_fonts()};assert all(fonts.values()),fonts
page.get_pixmap(matrix=fitz.Matrix(400/72,400/72),alpha=False).save(O/'figure1.png')
page.get_pixmap(matrix=fitz.Matrix(300/72,300/72),alpha=False).save(O/'preview.png')
ImageOps.grayscale(Image.open(O/'figure1.png')).save(O/'figure1_grayscale.png',dpi=(400,400))
with Image.open(O/'figure1.png') as im: im.copy().save(O/'figure1.png',dpi=(400,400))

# Bounded identity and export checks (same protected quantities as the source).
oldcoords=json.loads((O/'plot_coordinates.json').read_text())
assert oldcoords['all_20_pairs_shown'] and len(oldcoords['families'])==20
assert all(np.isfinite(v) for f in oldcoords['families'] for k in ['baseline_xy','copra_xy'] for v in f[k])
approved=R/'outputs/figure1_c_contours_20261031_v1/bootstrap_regions.json'
frozen_regions=json.loads((O/'bootstrap_regions.json').read_text())
def _assert_numclose(a,b,path='root',rtol=1e-9,atol=1e-12):
    """Recursively compare JSON values; tolerate last-bit BLAS/KDE round-off only."""
    if isinstance(a,dict):
        assert set(a)==set(b),'key mismatch at %s'%path
        for k in a:_assert_numclose(a[k],b[k],path+'/'+str(k),rtol,atol)
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b),'length mismatch at %s'%path
        for i,(x,y) in enumerate(zip(a,b)):_assert_numclose(x,y,path+'/%d'%i,rtol,atol)
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        np.testing.assert_allclose(a,b,rtol=rtol,atol=atol,err_msg=path)
    else:
        assert a==b,'value mismatch at %s: %r != %r'%(path,a,b)
if approved.exists():
    _assert_numclose(frozen_regions,json.loads(approved.read_text()),'bootstrap_regions')
base_pdf=R/'outputs/manuscript_layout_70pct_20261103/latex_source/figures/figure1.pdf'
text_check=None
if base_pdf.exists():
    from collections import Counter
    b=fitz.open(base_pdf);bn=Counter(w[4] for w in b[0].get_text('words'));b.close()
    nn=Counter(w[4] for w in page.get_text('words'))
    only_new=nn-bn;only_base=bn-nn
    text_check={'only_new':dict(list(only_new.items())[:20]),'only_baseline':dict(list(only_base.items())[:20]),
                'identical':not only_new and not only_base}
    assert text_check['identical'],text_check
summary={'status':'PASS','size_mm':[round(pw,2),round(ph,2)],'data_sha256':DATA_SHA256,
         'figure_text_identical_to_baseline':text_check,'embedded_fonts':fonts,
         'all_marks_vector':True,'raster_images':0,
         'n_pos_unseen_text_layer':{'chars':[hex(ord(c)) for c in added_chars],
             'native_sizes_pt':[7.5,5.25],
             'reason':'matplotlib >=3.10 draws non-BMP STIX math-italic glyphs as vector paths; '
                      'an invisible render-mode-3 text copy restores exact text extraction, '
                      'rendered pixels asserted unchanged'},
         'panel_a':{'columns_mm':[[COLS[0],COLS[0]+COLW],[COLS[1],COLS[1]+COLW],[COLS[2],COLS[2]+COLW]],
                    'scatter_y_mm':[20,49.5],'reduction_y_mm':[58,69],'native_geometry':True},
         'panel_b':{'matrix_x_mm':[MAT_X,MAT_X+MAT_W],'rows_mm':[MAT_TOP,MAT_TOP+56],
                    'column_width_mm':round(MAT_W/3,2),'row_height_mm':ROW_H,'n_pos_column_x_mm':N_POS_X},
         'panel_c':{'response_native_canvas_mm':[rp_w,rp_h],'crop_mm':[round(v,2) for v in clip_mm],
                    'crop_size_mm':[round(c_w,2),round(c_h,2)],'target_rect_mm':[round(v,2) for v in (c_tx0,c_ty0,c_tx0+c_w,c_ty0+c_h)],
                    'placement_scale':1.0,'fonts_preserved':True},
         'matrix_family_order':[r['id'] for r in rows],'matrix_column_scales':CMAX.tolist(),
         'matrix_mean':G.mean(axis=0).tolist(),'base_layout_audit':issues}
(O/'render_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
print(json.dumps({k:v for k,v in summary.items() if k not in ['embedded_fonts','matrix_family_order']},ensure_ascii=False))
