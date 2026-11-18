"""Three distinct PRA views: numerical map, beeswarm, paired radial SD profile."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,SymLogNorm,Normalize
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from PIL import Image
O=Path(__file__).resolve().parent;D=json.loads((O/'plot_data.json').read_text());P=json.loads((O/'restricted/pair_plot_data.json').read_text())
W,H=178,178
NAVY,LIGHT,BLUE,RED,INK,GREY,LINE='#3A3A3A','#F0F0F0','#56B4E9','#D55E00','#222222','#7A7A7A','#D5D5D5'
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.linewidth':.55,
    'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':2,'text.color':INK,'axes.labelcolor':INK,
    'xtick.color':INK,'ytick.color':INK,'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'dejavusans'})
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
PANEL='A'
def mapping(x,y,w=0,h=0):
    if PANEL=='A':return x*.71,y*1.36,w*.71,h*1.36
    if PANEL=='B':return 103+(x-3)*.91,y-121,w*.91,h
    if PANEL=='C':return 107+(x-85)*1.10,y-34,w*1.10,h
def axmm(x,y,w,h,**kw):
    x,y,w,h=mapping(x,y,w,h)
    return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H],**kw)
def text(x,y,s,**kw):
    x,y,_,_=mapping(x,y)
    return fig.text(x/W,1-y/H,s,va=kw.pop('va','top'),**kw)
def title(x,y,l,s):
    x,y,_,_=mapping(x,y)
    fig.text(x/W,1-y/H,l,fontsize=11,fontweight='bold',va='top')
    fig.text((x+5)/W,1-(y+.6)/H,s,fontsize=9,fontweight='bold',va='top')
def box(x,y,w,h,color,**kw):
    x,y,w,h=mapping(x,y,w,h)
    p=Rectangle((x/W,1-(y+h)/H),w/W,h/H,transform=fig.transFigure,facecolor=color,**kw);fig.add_artist(p)
def line(x1,y1,x2,y2,**kw):
    x1,y1,_,_=mapping(x1,y1)
    x2,y2,_,_=mapping(x2,y2)
    fig.add_artist(Line2D([x1/W,x2/W],[1-y1/H,1-y2/H],transform=fig.transFigure,**kw))
def style(ax):
    ax.spines[['top','right']].set_visible(False)
    for sp in ax.spines.values():sp.set_color(LINE)
    ax.tick_params(pad=2)
LABEL={'mcmc1':'mcmc1','NO1':'NO1','charged_float':'Charged residues','ele':'Electrostatics','no_count':'NO contacts','cc_count':'CC contacts',
    'torsional_clash':'Torsional clash','Incl1':'Inclination','Buckle1':'Buckle','rise1':'Rise'}
SHORT={'mcmc1':'mcmc','NO1':'NO1','charged_float':'Charge','ele':'Ele','no_count':'NO','cc_count':'CC','torsional_clash':'Clash','Incl1':'Incl','Buckle1':'Buckle','rise1':'Rise'}
allterms=[(b,t) for b,ss in D['branches'].items() for t in ss['terms']]
cmap=LinearSegmentedColormap.from_list('signed',['#00527F','#56B4E9','#FFFFFF','#E69F00','#B34700'])
norm=SymLogNorm(linthresh=.005,linscale=.65,vmin=-.7,vmax=.7,base=10)
assert all(abs(v)<=.7 for b in P.values() for r in b['rows'] for v in [*r['contributions'],r['delta_G']])

# A: complete exact numerical accounting with aligned display-state tracks.
title(3,4,'A','Term contributions and response')
text(8,13,'106 pairs · 121 units ·\n42 exact zeros · 13 additional rounding ties',fontsize=8,color=GREY)

def branch_map(branch,header,y0):
    pp=P[branch];ss=D['branches'][branch];rs=pp['rows'];n=len(rs);nf=len(pp['terms']);x=33;wid=99;cw=wid/n
    text(33,header,f'{branch} · {n} pairs',fontsize=8,fontweight='bold')
    text(33,header+3.8,f"{ss['exact_zero']}/{n} exact zeros · +{ss['rounding_only']} rounding "+('tie' if ss['rounding_only']==1 else 'ties'),fontsize=7.5,color=GREY)
    array=np.asarray([r['contributions'] for r in rs]).T
    ax=axmm(x,y0-1.35,wid,nf*3.2)
    ax.pcolormesh(np.arange(n+1),np.arange(nf+1),array,cmap=cmap,norm=norm,edgecolors='white',linewidth=.17)
    ax.set(xlim=(0,n),ylim=(nf,0));ax.axis('off')
    for i,k in enumerate(pp['terms']):
        y=y0+.25+3.2*i
        text(31,y,LABEL[k],ha='right',va='center',fontsize=7.5)
        text(134,y,f"{sum(r['changed'][i] for r in rs)}/{n}",va='center',fontsize=7.5,color=GREY)
    sy=y0+nf*3.2+2.2
    line(8,sy-2,140,sy-2,color=LINE,lw=.5)
    response=np.asarray([r['delta_G'] for r in rs]);ra=axmm(x,sy-1.3,wid,2.6)
    ra.pcolormesh(np.arange(n+1),[0,1],response[None,:],cmap=cmap,norm=norm,edgecolors='white',linewidth=.17);ra.axis('off')
    text(31,sy,r'$\sum_j\beta_j\Delta x_j$',fontsize=8,ha='right',va='center')
    for i,(label,key) in enumerate([('Exact zero','zero'),('2-decimal tie','rounded_zero'),('Added tie','rounding_only')]):
        y=sy+4.3+3.2*i;text(31,y,label,ha='right',va='center',fontsize=7.5)
        for j,r in enumerate(rs):
            col=LIGHT
            if r[key]:col=BLUE if r['rounding_only'] and i>0 else NAVY
            box(x+j*cw,y-1.2,cw*.95,2.4,col,lw=0)
        text(134,y,f"{sum(r[key] for r in rs)}/{n}",va='center',fontsize=7.5,color=BLUE if i==2 else GREY)
    nzero=sum(r['zero'] for r in rs)
    line(x+cw*nzero,y0-1.6,x+cw*nzero,y+1.3,color=BLUE,lw=.7)
    return ax
im=branch_map('ss+str',21,28);branch_map('ds+enz',64,71)
text(134,24,'n/N',fontsize=7.5,color=GREY)
text(134,67,'n/N',fontsize=7.5,color=GREY)
text(8,111,'Columns follow the same fixed feature-pattern order in every track.',fontsize=7.5,color=GREY)
ca=axmm(65,122,67,1.8);cb=fig.colorbar(im.collections[0],cax=ca,orientation='horizontal',ticks=[-.7,0,.7])
cb.outline.set_visible(False);cb.ax.tick_params(labelsize=7.5,pad=1,length=1.5);cb.ax.minorticks_off()
cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:g}'))
if cb.solids is not None:cb.solids.set_rasterized(False)
text(90,117,r'Term / response (kcal/mol)',ha='center',fontsize=7.5)

# B follows the visual grammar of the actual MPA-MutPred Fig. 4.
PANEL='B'
title(3,125,'B','Term contribution distributions')
text(8,134,'All pairs · ordered by mean absolute contribution',fontsize=7.5,color=GREY)
be=axmm(28,147,43,54);yy=[];labels=[]
wtmap=LinearSegmentedColormap.from_list('WTstep',['#DFF0F8','#9ED4EC','#56B4E9','#00527F'])
for bi,(branch,pp) in enumerate(P.items()):
    terms=pp['terms'];rs=pp['rows'];offset=0 if bi==0 else 4.8
    order=sorted(range(len(terms)),key=lambda j:-np.mean([abs(r['contributions'][j]) for r in rs]))
    for rank,j in enumerate(order):
        k=terms[j];y=offset+rank;yy.append(y);labels.append(LABEL[k].replace('Charged residues','Charged'))
        vals=np.asarray([r['contributions'][j] for r in rs]);colorvals=np.asarray([r['wt_percentiles'][j] for r in rs])
        bins=np.round(vals/.008).astype(int);shifts=np.zeros(len(vals))
        for tag in np.unique(bins):
            ids=np.flatnonzero(bins==tag);ids=np.random.default_rng(2014+j).permutation(ids)
            shifts[ids]=np.linspace(-.34,.34,len(ids)) if len(ids)>1 else 0
        be.scatter(vals,y+shifts,c=colorvals,cmap=wtmap,vmin=0,vmax=1,s=7,linewidths=0,alpha=.95,zorder=3)
        be.axhline(y,color='#DCDCDC',lw=.4,ls=(0,(1,2)),zorder=0)
be.axvline(0,color='#777777',lw=.65,zorder=1)
be.set(xlim=(-.42,.61),ylim=(10.5,-.7),xticks=[-.4,0,.4],yticks=yy,yticklabels=labels)
be.tick_params(axis='y',length=0,pad=3);style(be)
be.spines['left'].set_visible(False)
be.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:g}'))
be.set_xlabel('Term contribution (kcal/mol)',fontsize=7.5,labelpad=3)
text(8,141,'ss+str',fontsize=7.5,color=GREY)
text(8,168,'ds+enz',fontsize=7.5,color=GREY,bbox=dict(facecolor='white',edgecolor='none',pad=.2))
wa=axmm(75,150,1.2,48);wc=fig.colorbar(matplotlib.cm.ScalarMappable(norm=Normalize(0,1),cmap=wtmap),cax=wa,ticks=[0,1])
wc.ax.set_yticklabels(['Low','High']);wc.outline.set_visible(False);wc.ax.tick_params(labelsize=7.5,length=0,pad=1)
if wc.solids is not None:wc.solids.set_rasterized(False)
text(80,178,'WT feature value',rotation=90,fontsize=7.5,ha='center',va='center')

# C: merge the same four active-term pair ellipses, using their common physical units.
PANEL='C'
title(85,125,'C','Covariance geometry')
text(90,134,'Solid: between · dashed: within ×10',fontsize=7.5,color=GREY)
CV=json.loads((O/'covariance_data.json').read_text())
specs=[('ss+str','charged_float','ele','ss: Charged / Ele','#00527F',0.8),
       ('ds+enz','no_count','cc_count','ds: NO / CC','#56B4E9',1.0),
       ('ds+enz','no_count','torsional_clash','ds: NO / Clash','#B34700',1.2),
       ('ds+enz','cc_count','torsional_clash','ds: CC / Clash','#E69F00',1.4)]
theta=np.linspace(0,2*np.pi,361);unit=np.array([np.cos(theta),np.sin(theta)])
ax=axmm(97,140,44,34);report=[]
for b,k1,k2,label,col,lw in specs:
    c=CV[b];ix=[c['terms'].index(k1),c['terms'].index(k2)];entry=dict(branch=b,ordered_terms=[k1,k2],color=col)
    for mode,factor in [('between',1),('within',10)]:
        cov=np.asarray(c[mode])[np.ix_(ix,ix)]
        eig,vec=np.linalg.eigh(cov);assert eig.min()>-1e-12
        root=vec@np.diag(np.sqrt(np.maximum(eig,0)))
        np.testing.assert_allclose(root@root.T,cov,atol=1e-12,rtol=0)
        curve=(root@unit)*factor
        assert np.abs(curve[0]).max()<2.4 and np.abs(curve[1]).max()<4
        if mode=='between':ax.fill(curve[0],curve[1],color=col,alpha=.06,lw=0,zorder=1)
        ax.plot(curve[0],curve[1],color=col,lw=lw if mode=='between' else lw+.25,
                ls='-' if mode=='between' else (0,(3,1.7)),alpha=.9,zorder=2 if mode=='between' else 4)
        entry[mode]=cov.tolist()
    report.append(entry)
ax.axhline(0,color=LINE,lw=.6,zorder=0);ax.axvline(0,color=LINE,lw=.6,zorder=0)
ax.set(xlim=(-2.4,2.4),ylim=(-4,4),xticks=[-2,0,2],yticks=[-3,0,3])
ax.set_xlabel('First term (kcal/mol)',fontsize=7.5,labelpad=2)
ax.set_ylabel('Second term (kcal/mol)',fontsize=7.5,labelpad=2)
style(ax);ax.tick_params(labelsize=7.5,pad=2)
fig.legend(handles=[Line2D([],[],color=GREY,lw=.9,label='Between'),Line2D([],[],color=GREY,lw=1.1,ls='--',label='Within ×10')],
    loc='upper center',bbox_to_anchor=(139/W,1-150/H),ncol=2,frameon=False,fontsize=7.5,
    borderaxespad=0,handlelength=1.3,handletextpad=.35,columnspacing=.8)
fig.legend(handles=[Line2D([],[],color=s[4],lw=s[5],label=s[3]) for s in specs],
    loc='upper center',bbox_to_anchor=(139/W,1-157/H),ncol=1,frameon=False,fontsize=7.5,
    borderaxespad=0,handlelength=1.5,handletextpad=.5,labelspacing=.3)
fig.text(6/W,1-172/H,'C: colors identify ordered term pairs; within ellipses use ×10 linear magnification.',fontsize=7.5,color=GREY,va='top')
(O/'ellipse_data.json').write_text(json.dumps(dict(physical_units='kcal/mol',coordinate_normalization='none',
    within_linear_display_multiplier=10,all_four_original_active_pairs_retained=True,pairs=report),indent=2))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from visual_qa import audit_layout
issues=audit_layout(fig)
(O/'layout_audit.json').write_text(json.dumps(dict(issues=issues,size_mm=[W,H]),indent=2))
fig.savefig(O/'figure4.png',dpi=400,facecolor='white')
for ext in ['pdf','svg']:fig.savefig(O/f'figure4.{ext}',facecolor='white')
Image.open(O/'figure4.png').convert('L').save(O/'figure4_grayscale.png',dpi=(400,400))
print(json.dumps(dict(outputs=['figure4.png','figure4.pdf','figure4.svg'],layout_issues=issues,size_mm=[W,H])))
