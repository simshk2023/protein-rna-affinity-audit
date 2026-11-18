"""Amplitude/rank annular plot for the 178 mm wide Figure 1 reflow.

Native point sizes are preserved (7.5/8/8.5 pt); only the canvas and the axes
rectangle grow from the original 82.4 x 76 mm to 104 x 88 mm so that panel c
can fill the wider lower-right region. Data, annular transform, bootstrap
regions and every label are identical to
figure1_bright_semicircle_20261101_v3/render_response.py.
"""
from pathlib import Path
import json,sys,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from PIL import Image,ImageOps
from scipy.stats import gaussian_kde

O=Path(__file__).resolve().parent
DATA_SHA256='6dfeefd6cc8a1469fc227103ad2713b7135e7cd8d18bf02f5a5eba66dd866e8a'
assert hashlib.sha256((O/'family_aggregates.json').read_bytes()).hexdigest()==DATA_SHA256,'Frozen family_aggregates.json changed.'
d=json.loads((O/'family_aggregates.json').read_text());rows=d['local']
from bright_palette import BLUE, TEAL, CORAL, AMBER, INK, SECONDARY, AXIS, GRID, LIGHT_BLUE, LIGHT_TEAL
GRAY=SECONDARY
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'text.color':INK,'pdf.fonttype':42,'svg.fonttype':'none'})
W,H=104,88
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
ax=fig.add_axes([.13,.34,.74,.58]);ax.set_aspect('equal');ax.axis('off')
ax.set_xlim(-2.92,2.92);ax.set_ylim(-.48,3.12)

QMIN=.005;QMAX=5.;BASE=.6;FACTOR=.65
def radius(q):return BASE+FACTOR*np.log10(np.asarray(q)/QMIN)
def project(q,rho):
    rr=radius(q);th=np.arccos(np.asarray(rho))
    return np.stack([rr*np.cos(th),rr*np.sin(th)],axis=-1)
def curve(q1,rho1,q2,rho2,n=120):
    # Interpolate only the drawing path; there are no intermediate model scores.
    t=np.linspace(0,1,n);q=10**((1-t)*np.log10(q1)+t*np.log10(q2))
    th=(1-t)*np.arccos(rho1)+t*np.arccos(rho2)
    rr=radius(q);return np.column_stack([rr*np.cos(th),rr*np.sin(th)])
THETA_MAX=np.pi
theta=np.linspace(0,THETA_MAX,361);RMAX=float(radius(QMAX))

# Sparse guide grid. Ratio 1 is highlighted; no RMS contours or ideal-error target.
for q in [.01,.1,1,3]:
    rr=float(radius(q));ax.plot(rr*np.cos(theta),rr*np.sin(theta),color='#51BAE7' if q==1 else GRID,lw=.95 if q==1 else .45,ls=(0,(4,3)) if q==1 else (0,(2,3)),zorder=0)
    ax.plot([rr,rr],[0,-.045],color=INK,lw=.55)
    if q not in [.03,.3]:ax.text(rr,-.11,f'{q:g}',fontsize=7.5,ha='center',va='top')
for rr in [BASE,RMAX]:
    ax.plot(rr*np.cos(theta),rr*np.sin(theta),color=AXIS,lw=.7,zorder=0)
ax.plot([-RMAX,RMAX],[0,0],color=AXIS,lw=.7,zorder=0)
for rho in [-1,-.5,0,.2,.4,.6,.8,1]:
    th=np.arccos(rho);co=np.cos(th);si=np.sin(th)
    ax.plot([BASE*co,RMAX*co],[BASE*si,RMAX*si],color=GRID,lw=.45,ls=(0,(2,3)),zorder=0)
    ax.plot([RMAX*co,(RMAX+.035)*co],[RMAX*si,(RMAX+.035)*si],color=INK,lw=.6)
    ax.text((RMAX+.135)*co,(RMAX+.135)*si,f'{rho:g}'.replace('-','−'),fontsize=7.5,ha='center',va='center')
ax.text(0,2.99,r'Spearman $\rho$',fontsize=8.5,ha='center',va='center')
ax.text(1.56,-.72,'SD ratio (log radius)',fontsize=8,ha='center',va='center')

# Joint family resampling in original metric space, before display remapping.
# Curves are bootstrap mean regions, NOT Taylor centred-error contours.
NBOOT=10000;SEED=20260921
draw=np.random.default_rng(SEED).integers(0,len(rows),(NBOOT,len(rows)))
boot={m:np.array([r[m][3:] for r in rows])[draw].mean(axis=1) for m in ['baseline','copra']}
region_audit={}
for m,style in [('baseline','--'),('copra','-')]:
    sample=boot[m];kde=gaussian_kde(sample.T)
    sd=np.sqrt(np.diag(kde.covariance))
    low=sample.min(axis=0)-4*sd;high=sample.max(axis=0)+4*sd
    low=np.maximum(low,[QMIN+.00001,-.599]);high=np.minimum(high,[QMAX-.00001,.999])
    qs=np.linspace(low[0],high[0],180);rs=np.linspace(low[1],high[1],180)
    Q,P=np.meshgrid(qs,rs);Z=kde(np.vstack([Q.ravel(),P.ravel()])).reshape(Q.shape)
    cell=(qs[1]-qs[0])*(rs[1]-rs[0]);mass=Z.sum()*cell
    assert .99<mass<1.01,(m,mass)
    ordered=np.sort(Z.ravel())[::-1];cumulative=np.cumsum(ordered)/ordered.sum()
    thresholds={str(p):float(ordered[np.searchsorted(cumulative,p)]) for p in [.5,.95]}
    xy=project(Q,P)
    ax.contour(xy[:,:,0],xy[:,:,1],Z,levels=[thresholds['0.95'],thresholds['0.5']],
               colors=[TEAL if m=='baseline' else BLUE]*2,linestyles=style,linewidths=[.85,.55],alpha=.75,zorder=2)
    region_audit[m]={'estimated_density_grid_mass':float(mass),'thresholds':thresholds,
                     'probability_mass_above_threshold':{p:float(Z[Z>=v].sum()/Z.sum()) for p,v in thresholds.items()},
                     'original_metric_bounds':[low.tolist(),high.tolist()]}

export=[]
for row in rows:
    aq,ar=row['baseline'][3:];cq,cr=row['copra'][3:]
    assert QMIN<aq<QMAX and QMIN<cq<QMAX and -.6<ar<1 and -.6<cr<1
    path=curve(aq,ar,cq,cr)
    ax.plot(path[:,0],path[:,1],color=LIGHT_BLUE,lw=.5,alpha=.55,zorder=1)
    # The head identifies the destination; markers encode the endpoint models.
    ax.add_patch(FancyArrowPatch(path[-9],path[-1],arrowstyle='-|>',mutation_scale=4.8,lw=.55,color=LIGHT_BLUE,alpha=.75,zorder=1))
    a=project(aq,ar);c=project(cq,cr)
    ax.scatter(*a,s=11,marker='o',facecolor='white',edgecolor=TEAL,lw=.6,zorder=3)
    ax.scatter(*c,s=15,marker='^',facecolor=BLUE,edgecolor='white',lw=.25,zorder=4)
    export.append({'family':row['id'],'baseline_sd_ratio':aq,'baseline_spearman':ar,'copra_sd_ratio':cq,'copra_spearman':cr,'baseline_xy':a.tolist(),'copra_xy':c.tolist()})

# Macro means are computed on the original score scale, before nonlinear mapping.
mean={m:np.mean([r[m][3:] for r in rows],axis=0) for m in ['baseline','copra']}
path=curve(*mean['baseline'],*mean['copra'])
ax.plot(path[:,0],path[:,1],color=AMBER,lw=1.35,zorder=5)
ax.add_patch(FancyArrowPatch(path[-15],path[-1],arrowstyle='-|>',mutation_scale=7,lw=1,color=AMBER,zorder=5))
ax.scatter(*project(*mean['baseline']),s=26,marker='o',facecolor='white',edgecolor=AMBER,lw=1.0,zorder=6)
ax.scatter(*project(*mean['copra']),s=32,marker='^',facecolor=AMBER,edgecolor='white',lw=.4,zorder=7)
# Reference-like right key: model symbols, contour styles and reference guide.
key=fig.add_axes([.04,.10,.92,.19]);key.axis('off')
key.legend(handles=[
    Line2D([],[],color=TEAL,ls='--',lw=.75,marker='o',mfc='white',mec=TEAL,ms=4,label='Additive'),
    Line2D([],[],color=BLUE,ls='-',lw=.75,marker='^',mfc=BLUE,mec='white',mew=.3,ms=4.8,label='CoPRA'),
    Line2D([],[],color=AMBER,lw=1.3,label='Mean shift'),
    Line2D([],[],color='#51BAE7',ls=(0,(4,3)),lw=.95,label='SD ratio = 1')],
    loc='upper center',ncol=2,frameon=False,fontsize=7.5,handletextpad=.55,labelspacing=1.0,columnspacing=1.6,handlelength=1.5,borderaxespad=0)
fig.text(.5,.105,'Contours: 50% / 95% bootstrap mean regions',fontsize=7.5,ha='center',va='bottom',color=GRAY)

change=boot['copra']-boot['baseline']
regions={'resamples':NBOOT,'seed':SEED,
    'unit':'whole family, same indices for both models and both metrics',
    'region_definition':'50/95% highest-density sets of Gaussian KDE of bootstrap means in original SD-ratio/Spearman coordinates; then mapped to annular coordinates',
    'kde_bandwidth':'Scott default','grid_shape':[180,180],'region_audit':region_audit,
    'paired_gain_estimate':(mean['copra']-mean['baseline']).tolist(),
    'bootstrap_gain_distribution_mean':change.mean(axis=0).tolist(),
    'paired_gain_ci95':np.quantile(change,[.025,.975],axis=0).tolist()}

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

# The approved bootstrap_regions.json is frozen: never overwrite it, but prove
# that the re-rendered contours match it to floating-point round-off.
FROZEN=O/'bootstrap_regions.json'
if FROZEN.exists():
    _assert_numclose(regions,json.loads(FROZEN.read_text()),'bootstrap_regions')
else:
    FROZEN.write_text(json.dumps(regions,indent=2))

# Inverse-coordinate checks protect against swapped metrics or an altered transform.
for row in export:
    for m in ['baseline','copra']:
        x,y=row[m+'_xy'];rr=np.hypot(x,y)
        np.testing.assert_allclose([QMIN*10**((rr-BASE)/FACTOR),x/rr],[row[m+'_sd_ratio'],row[m+'_spearman']],atol=1e-12,rtol=1e-12)
sys.path.insert(0,str(Path(__file__).resolve().parent))
from visual_qa import audit_layout
issues=audit_layout(fig);(O/'response_wide_layout_audit.json').write_text(json.dumps(issues,indent=2))
assert not issues,issues
fig.savefig(O/'response_wide.png',dpi=450,facecolor='white')
ImageOps.grayscale(Image.open(O/'response_wide.png')).save(O/'response_wide_grayscale.png')
for ext in ['pdf','svg']:fig.savefig(O/f'response_wide.{ext}',facecolor='white')
(O/'plot_coordinates.json').write_text(json.dumps({'family_aggregates_only':True,'families':export,'mean_scores_before_mapping':{m:a.tolist() for m,a in mean.items()},'mapping':{'display_radius':'0.6 + 0.65*log10(SD_ratio/0.005)','display_angle':'arccos(position_Spearman)','amplitude_bounds':[QMIN,QMAX],'rank_bounds':[-1,1]},'size_mm':[W,H],'frame':'full semicircle','angle_ticks':[-1,-.5,0,.2,.4,.6,.8,1],'labelled_radius_ticks':[.01,.1,1,3],'data_source_sha256':DATA_SHA256,'all_20_pairs_shown':len(export)==20,'inverse_mapping_check':'PASS','layout_audit':issues},indent=2))
print(json.dumps({'families':len(rows),'size_mm':[W,H],'frame':'full semicircle','angle_ticks':[-1,-.5,0,.2,.4,.6,.8,1],'labelled_radius_ticks':[.01,.1,1,3],'bootstrap_region_audit':region_audit,'layout_audit':issues}))
