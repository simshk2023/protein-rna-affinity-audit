"""Matched mutation-effect MSE decomposition: main, uncertainty, accounting, intervention."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Rectangle,ConnectionPatch,FancyBboxPatch,FancyArrowPatch,Circle
from matplotlib.colors import SymLogNorm,LinearSegmentedColormap
from scipy.stats import gaussian_kde
from PIL import Image

O=Path(__file__).resolve().parent;D=json.loads((O/'plot_data.json').read_text())
W,H=178,188
BLUE,SEED,RED,INK,AXIS,GRID='#4DBBD5','#91D1C2','#E64B35','#252525','#777777','#E3E9EC'
TEAL,INDIGO='#00A087','#3C5488'
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,
    'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.linewidth':.6,
    'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2.3,'ytick.major.size':2.3,
    'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'dejavusans',
    'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
PANEL='A'
def mapping(x,y,w=0,h=0):
    if PANEL=='A':return x,y,w,h
    if PANEL=='B':return 3+(x-83)*1.15,y+93,w*1.15,h
    if PANEL=='C':return 93+(x-3)*.59,y-92,w*.59,h
    if PANEL=='DG':return x+90,y-93,w,h
    if PANEL=='DO':return 109+(x-88)*1.15,y-43,w*1.15,h
def rect(x,y,w,h):
    x,y,w,h=mapping(x,y,w,h)
    return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def txt(x,y,s,**kw):
    x,y,_,_=mapping(x,y)
    return fig.text(x/W,1-y/H,s,va=kw.pop('va','top'),**kw)
def title(x,y,l,s):
    x,y,_,_=mapping(x,y)
    fig.text(x/W,1-y/H,l,fontsize=11,fontweight='bold',va='top')
    fig.text((x+5)/W,1-(y+.6)/H,s,fontsize=9,fontweight='bold',va='top')
def style(ax,box=False):
    for spine in ax.spines.values():spine.set_color(AXIS);spine.set_linewidth(.55)
    if not box:ax.spines[['top','right']].set_visible(False)
    ax.tick_params(pad=2)
def zero(ax):
    ax.axvline(0,color=AXIS,lw=.55,zorder=0);ax.axhline(0,color=AXIS,lw=.55,zorder=0)
g=np.array([f['gain'] for f in D['main']]);s=np.array([f['seed_gains'] for f in D['main']]);mean=g.mean(axis=0)
def scatter_a(ax,small=False):
    zero(ax);ax.plot([-.1,1.1],[.1,-1.1],color=AXIS,lw=.7,ls=(0,(3,2)),zorder=0)
    ax.scatter(s[:,:,0].ravel(),s[:,:,1].ravel(),s=9 if small else 13,
               facecolors='none',edgecolors=SEED,linewidths=.6,alpha=.9,zorder=2)
    ax.scatter(g[:,0],g[:,1],s=17 if small else 24,color=BLUE,edgecolors='white',linewidths=.4,zorder=3)
    ax.scatter(mean[0],mean[1],s=32 if small else 45,marker='D',color=RED,edgecolors='white',linewidths=.5,zorder=5)

# A: primary family-level analysis, with a mapped and explicitly bounded zoom.
title(3,4,'A','Mean–local decomposition')
txt(16,12,'220 pairs · 20 families · 3 seeds',fontsize=7.5,color=AXIS)
legend=[Line2D([],[],ls='',marker='o',mfc=BLUE,mec='white',ms=4,label='Family'),
        Line2D([],[],ls='',marker='o',mfc='none',mec=SEED,ms=3.5,label='Seed'),
        Line2D([],[],ls='',marker='D',mfc=RED,mec='white',ms=4.5,label='Overall')]
fig.legend(handles=legend,loc='upper left',bbox_to_anchor=(16/W,1-18/H),ncol=3,
           fontsize=7.5,frameon=False,borderaxespad=0,handletextpad=.3,columnspacing=.9)
a=rect(16,25,59,44);scatter_a(a)
a.set(xlim=(-.08,.9),ylim=(-1.05,.18),xticks=[0,.3,.6,.9],yticks=[-1,-.5,0])
a.set_xlabel(r'Mean-component reduction (pKd$^2$)',labelpad=3)
a.set_ylabel(r'Local-component reduction (pKd$^2$)',labelpad=3);style(a)
worst=int(np.argmin(g[:,1]));a.text(g[worst,0]+.025,g[worst,1]-.015,D['main'][worst]['family'],fontsize=7.5,color=BLUE)
box=(-.02,.23,-.055,.12)
a.add_patch(Rectangle((box[0],box[2]),box[1]-box[0],box[3]-box[2],fill=False,ec=BLUE,lw=.6,ls=(0,(2,2))))
zoom=rect(36,43,25,19);scatter_a(zoom,small=True)
zoom.set(xlim=box[:2],ylim=box[2:],xticks=[0,.1,.2],yticks=[0,.1]);style(zoom,box=True)
zoom.set_title('Near-zero region',fontsize=7.5,pad=3,
               bbox=dict(facecolor='white',edgecolor='none',pad=.3))
fig.add_artist(ConnectionPatch(xyA=(box[1],box[3]),coordsA=a.transData,xyB=(0,1),coordsB=zoom.transAxes,color=GRID,lw=.6,zorder=1))
txt(16,86,r'$\mathrm{MSE}_{\Delta}=\mu^2+\langle(q-\mu)^2\rangle$',fontsize=8)

# B: the joint bootstrap distribution preserves paired component uncertainty.
PANEL='B'
title(83,4,'B','Joint uncertainty')
txt(88,14,'10,000 paired family draws',fontsize=7.5,color=AXIS)
jb=D['joint_bootstrap'];bx=np.asarray(jb['x']);by=np.asarray(jb['y']);z=np.asarray(jb['density'])
ordered=np.sort(z.ravel())[::-1];cs=np.cumsum(ordered)/ordered.sum()
levels=[ordered[np.searchsorted(cs,p)] for p in [.95,.8,.5]]+[z.max()*1.0001]
b=rect(98,25,35,34)
b.contourf(bx,by,z,levels=levels,colors=['#E5F5F9','#A6DDEA','#4DBBD5'])
b.contour(bx,by,z,levels=levels[:-1],colors=BLUE,linewidths=.55)
zero(b);b.plot([-.04,.23],[.04,-.23],color=AXIS,lw=.75,ls=(0,(3,2)))
b.scatter(mean[0],mean[1],marker='D',s=33,color=RED,edgecolors='white',linewidths=.5,zorder=4)
ci_x=D['effects']['mean_mse']['ci95'];ci_y=D['effects']['centered_mse']['ci95']
b.plot(ci_x,[mean[1]]*2,color=BLUE,lw=.7,zorder=3)
b.plot([mean[0]]*2,ci_y,color=BLUE,lw=.7,zorder=3)
b.set(xlim=(bx[0],bx[-1]),ylim=(by[0],by[-1]),xticks=[0,.1,.2,.3],yticks=[-.2,-.1,0])
b.set_xlabel(r'Mean reduction (pKd$^2$)',fontsize=7.5,labelpad=2)
b.set_ylabel(r'Local reduction (pKd$^2$)',fontsize=7.5,labelpad=2);style(b)
mt=rect(98,19,35,4);mt.fill_between(bx,z.sum(axis=0),color=SEED,lw=0);mt.set_xlim(bx[0],bx[-1]);mt.axis('off')
mr=rect(134.5,25,3.5,34);mr.fill_betweenx(by,z.sum(axis=1),color=RED,alpha=.4,lw=0);mr.set_ylim(by[0],by[-1]);mr.axis('off')
txt(87,74,'Estimate [95% CI]',fontsize=7.5,color=AXIS)
for yy,key,label in [(78.5,'mean_mse','Mean'),(83,'centered_mse','Local'),(87.5,'total_mse','Total')]:
    v=D['effects'][key];lo,hi=v['ci95']
    txt(87,yy,label,fontsize=7.5)
    txt(99,yy,f"{v['estimate']:+.3f} [{lo:+.3f}, {hi:+.3f}]",fontsize=7.5)

# C: actual family-by-seed additive costs on a common signed scale.
PANEL='C'
title(3,96,'C','Family error accounting')
txt(24,104,'Families ordered by local-MSE increase',fontsize=7.5,color=AXIS)
txt(24,108,'Positive = more error',fontsize=7.5,color=RED)
atlas=D['accounting_atlas'];order=np.asarray(atlas['family_order']);matrix=np.asarray(atlas['matrix'])
c=rect(24,117,100,18)
norm=SymLogNorm(linthresh=.01,linscale=.6,vmin=-1,vmax=1,base=10)
cmap=LinearSegmentedColormap.from_list('cost',[BLUE,'#B8E5EF','#FFFFFF','#F39B7F',RED])
im=c.pcolormesh(np.arange(21)-.5,np.arange(10)-.5,matrix,cmap=cmap,norm=norm,linewidth=.2,edgecolors='white')
c.set(xlim=(-.5,19.5),ylim=(8.5,-.5),xticks=np.arange(20),xticklabels=[D['main'][i]['family'] for i in order],yticks=[])
c.xaxis.tick_top();c.tick_params(axis='x',length=0,pad=2,labelsize=7.5)
plt.setp(c.get_xticklabels(),rotation=90,ha='center',va='bottom')
for spine in c.spines.values():spine.set_visible(False)
for y in [2.5,5.5]:c.axhline(y,color=AXIS,lw=.6)
for yy,label in [(120,r'$\Delta V_{\mathrm{pred}}$'),(126,r'$-2\Delta C$'),(132,r'$\Delta\mathrm{MSE}_{\mathrm{local}}$')]:
    txt(21,yy,label,ha='right',va='center',fontsize=8)
txt(133,112,'Mean',ha='center',fontsize=7.5,fontweight='bold')
summary=[np.mean(matrix[:3]),np.mean(matrix[3:6]),np.mean(matrix[6:])]
for yy,value in zip([120,126,132],summary):txt(133,yy,f'{value:+.4f}',ha='center',va='center',fontsize=7.5,color=RED if value>0 else BLUE)
txt(24,140,'Three strips per family: S1 / S2 / S3',fontsize=7.5)
ca=rect(100,141,24,1.5);cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-1,0,1])
cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda v,pos:f'{v:g}'));cb.ax.minorticks_off()
cb.ax.tick_params(labelsize=7.5,length=1.5,pad=1);cb.outline.set_visible(False)
if cb.solids is not None:cb.solids.set_rasterized(False)
txt(112,137,r'Change (pKd$^2$)',fontsize=7.5,ha='center')

# D: two integrated views; matched residual geometry and 2D outcome trajectories.
PANEL='DG'
title(3,151,'D','Supervision geometry and validation')
txt(16,159,'Auxiliary losses',fontsize=8,fontweight='bold')

txt(16,165,r'$L_B=\frac{1}{2}(e_W^2+e_M^2)$',fontsize=8,color=INDIGO)
txt(49,165,r'$L_C=\frac{1}{2}(e_M-e_W)^2$',fontsize=8,color=RED)
u=np.linspace(-1.5,1.5,241);xx,yy=np.meshgrid(u,u)
lossB=.5*(xx**2+yy**2);lossC=.5*(yy-xx)**2
geo=rect(20,175,40,32)
# Exact auxiliary-loss contrast: C minus B = -eW*eM.
contrast=lossC-lossB
np.testing.assert_allclose(contrast,-xx*yy,atol=1e-14,rtol=0)
bgmap=LinearSegmentedColormap.from_list('contrast',[BLUE,'#D6EFF5','white','#FADED8',RED])
bg=geo.contourf(xx,yy,contrast,levels=np.linspace(-2.25,2.25,25),cmap=bgmap,alpha=.42,vmin=-2.25,vmax=2.25,zorder=0)
levels=[.0625,.125,.25,.5,1,2]
bc=geo.contour(xx,yy,lossB,levels=levels,colors=BLUE,linewidths=.75,zorder=2)
cc=geo.contour(xx,yy,lossC,levels=levels,colors=RED,linewidths=.65,linestyles='dashed',zorder=2)
geo.clabel(bc,levels=[.25,1],fmt='%g',fontsize=7.5,inline=True)
geo.plot(u,u,color=RED,lw=1.1,alpha=.65,zorder=1)
geo.scatter([0],[0],color=BLUE,edgecolor=INDIGO,s=16,zorder=5,linewidth=.5)
# Matched analytic descent steps, not observed parameter-training trajectories.
anchors=np.array([[-1.05,.4],[.7,1.05],[-.65,-1.05],[1.05,-.55]])
eta=.24
for px,py in anchors:
    point=np.array([px,py]); gradB=point; gradC=np.array([px-py,py-px])
    geo.scatter(px,py,s=7,facecolors='white',edgecolors=INDIGO,linewidths=.4,zorder=5)
    for grad,col,ls in [(gradB,BLUE,'-'),(gradC,RED,'--')]:
        end=point-eta*grad
        assert np.max(np.abs(end))<1.5
        old=.5*(px*px+py*py) if col==BLUE else .5*(py-px)**2
        new=.5*(end@end) if col==BLUE else .5*(end[1]-end[0])**2
        assert new<old
        geo.annotate('',xy=end,xytext=point,arrowprops=dict(arrowstyle='-|>',color='white',lw=2,shrinkA=0,shrinkB=0,mutation_scale=6),zorder=6)
        geo.annotate('',xy=end,xytext=point,arrowprops=dict(arrowstyle='-|>',color=col,lw=.85,ls=ls,shrinkA=0,shrinkB=0,mutation_scale=6),zorder=7)
geo.set(xlim=(-1.5,1.5),ylim=(-1.5,1.5),xticks=[-1,0,1],yticks=[-1,0,1],aspect='equal')
geo.set_xlabel(r'$e_W$ (pKd)',fontsize=7.5,labelpad=1)
geo.set_ylabel(r'$e_M$ (pKd)',fontsize=7.5,labelpad=1);style(geo,box=True)
geo.axhline(0,color=AXIS,lw=.35,alpha=.6,zorder=1);geo.axvline(0,color=AXIS,lw=.35,alpha=.6,zorder=1)
colorax=rect(64,178,1.4,25)
cb=fig.colorbar(bg,cax=colorax,ticks=[-2,0,2]);cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=7.5,length=1.5,pad=1)
if cb.solids is not None:cb.solids.set_rasterized(False)
txt(65,173,r'$L_C-L_B$',fontsize=7.5,ha='center')
fig.legend(handles=[Line2D([],[],color=BLUE,lw=1,label='Endpoint B'),Line2D([],[],color=RED,lw=1,ls='--',label='Difference C')],
    loc='upper center',bbox_to_anchor=(130/W,1-121/H),ncol=2,frameon=False,fontsize=7.5,borderaxespad=0,
    handlelength=1.2,handletextpad=.3,columnspacing=.7)

PANEL='DO'
txt(88,170,'54 pairs · 4 families · 3 seeds',fontsize=7.5,color=AXIS)
armB=np.asarray(D['supervision_arm_components']['B']);armC=np.asarray(D['supervision_arm_components']['C'])
vg=armB-armC;vm=vg.mean(axis=1);familycols=[BLUE,TEAL,RED,INDIGO];familymarks=['o','s','^','D']
out=rect(88,175,46,32)
for i,(col,mark) in enumerate(zip(familycols,familymarks)):
    for si in range(3):
        p=armB[i,si,:2];q=armC[i,si,:2]
        out.plot([p[0],q[0]],[p[1],q[1]],color=col,alpha=.35,lw=.65,zorder=1)
        out.scatter(*p,s=10,marker='o',facecolors='white',edgecolors=col,linewidths=.6,zorder=2)
        out.scatter(*q,s=10,marker='s',facecolors=col,edgecolors='white',linewidths=.3,zorder=3)
    p=armB[i,:,:2].mean(axis=0);q=armC[i,:,:2].mean(axis=0)
    out.annotate('',xy=q,xytext=p,arrowprops=dict(arrowstyle='->',color=col,lw=1,shrinkA=0,shrinkB=0,mutation_scale=5),zorder=4)
    out.scatter(*p,s=23,marker='o',facecolors='white',edgecolors=col,linewidths=.8,zorder=5)
    out.scatter(*q,s=21,marker='s',facecolors=col,edgecolors='white',linewidths=.4,zorder=6)
    out.annotate(f'V{i+1:02d}',xy=p,xytext=([5,6] if i!=2 else [5,-9]),textcoords='offset points',fontsize=7.5,color=col)
out.set(xlim=(-.1,3.15),ylim=(.125,.35),xticks=[0,1,2,3],yticks=[.15,.25,.35])
out.set_xlabel(r'Mean MSE (pKd$^2$)',fontsize=7.5,labelpad=1)
out.set_ylabel(r'Local MSE (pKd$^2$)',fontsize=7.5,labelpad=1);style(out)
out.grid(color=GRID,lw=.4,zorder=0)
fig.text(159/W,1-127/H,'○ B → ■ C',fontsize=7.5,color=INK,va='top')

# The inset displays paired reductions, with its own explicit axes and units.
delta=rect(95,178,21,12)
for i,(col,mark) in enumerate(zip(familycols,familymarks)):
    delta.scatter(vg[i,:,0],vg[i,:,1],s=10,marker=mark,facecolors='white',edgecolors=col,linewidths=.65,zorder=3)
    delta.scatter(vm[i,0],vm[i,1],s=20,marker=mark,c=col,edgecolors='white',linewidths=.3,zorder=4)
zero(delta);delta.set(xlim=(-.08,.13),ylim=(-.018,.012),xticks=[-.05,.1],yticks=[-.01,.01])
delta.tick_params(labelsize=7.5,pad=1,length=1.5)
delta.set_title('B − C reductions',fontsize=7.5,pad=2,
    bbox=dict(facecolor='white',edgecolor='none',pad=.1))
delta.set_xlabel('Mean',fontsize=7.5,labelpad=0)
delta.set_ylabel('Local',fontsize=7.5,labelpad=0);style(delta,box=True)
txt(88,215,r'Both retain absolute-affinity loss',fontsize=7.5,color=AXIS)

fig.canvas.draw();renderer=fig.canvas.get_renderer()
pts=np.vstack([g[:,:2],s[:,:,:2].reshape(-1,2)])
covered=[i for i,p in enumerate(a.transData.transform(pts)) if zoom.get_window_extent(renderer).contains(*p)]
assert not covered
np.testing.assert_allclose(vg,np.asarray([r['seed_gains'] for r in D['supervision']['families']]),atol=1e-12,rtol=0)
allout=np.concatenate([armB[:,:,:2].reshape(-1,2),armC[:,:,:2].reshape(-1,2)])
outcovered=[i for i,p in enumerate(out.transData.transform(allout)) if delta.get_window_extent(renderer).contains(*p)]
assert not outcovered,outcovered
assert np.all((allout[:,0]>out.get_xlim()[0])&(allout[:,0]<out.get_xlim()[1])&(allout[:,1]>out.get_ylim()[0])&(allout[:,1]<out.get_ylim()[1]))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from visual_qa import audit_layout
issues=audit_layout(fig)
(O/'layout_audit.json').write_text(json.dumps(dict(issues=issues,size_mm=[W,H],all_points_inside_axes=True,
    inset_covers_no_points=True,D_inset_covers_no_points=True),indent=2))
fig.savefig(O/'figure3.png',dpi=400,facecolor='white')
for ext in ['pdf','svg']:fig.savefig(O/f'figure3.{ext}',facecolor='white')
Image.open(O/'figure3.png').convert('L').save(O/'figure3_grayscale.png',dpi=(400,400))
print(json.dumps(dict(outputs=['figure3.png','figure3.pdf','figure3.svg'],layout_issues=issues,size_mm=[W,H])))
