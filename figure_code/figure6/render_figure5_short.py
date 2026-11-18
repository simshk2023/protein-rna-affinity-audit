"""Wide-format Figure 5 (178 x 158 mm): same panels/data/typography, recomposed for the 185.6 mm text column."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch,Rectangle
from matplotlib.ticker import NullLocator
from scipy.stats import gaussian_kde
from PIL import Image
O=Path(__file__).resolve().parent;d=json.loads((O/'plot_data.json').read_text());dd=json.loads((O/'plot_data_D.json').read_text());v=dd['values']
G='#56B4E9';DARK='#2B2B2B';P='#D55E00';PD='#8A3A00';GRAY='#7A7A7A';INK='#222222';CB1='#0072B2';CB2='#E69F00'
plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':'#555859','xtick.color':INK,'ytick.color':INK,'axes.linewidth':.55,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','xtick.major.width':.5,'ytick.major.width':.5,'hatch.linewidth':.55})
# Short wide canvas: 178 x 143 mm (for two figures on one page at 70% display). Font sizes unchanged.
f=plt.figure(figsize=(178/25.4,135/25.4),facecolor='white')
def txt(x,y,s,**kw):return f.text(x,y,s,**kw)
def clean(a):
 a.spines[['top','right']].set_visible(False);a.tick_params(length=2.5)
for x,letter in [(.025,'A'),(.548,'B')]:txt(x,.969,letter,fontsize=11,weight='bold')
# A: accepted checkpoint matrix, with aligned edges and contrast-aware text.
a=f.add_axes([.105,.586,.395,.351]);models=['0','1','2','3','4','ensemble'];mat=np.array([d['checkpoints'][k]['ratios']+[d['checkpoints'][k]['mean']] for k in models]).T
cm=LinearSegmentedColormap.from_list('response',['#EAF5FB','#BFE3F3','#7FC6EA','#00527F']);norm=LogNorm(1e-5,1)
im=a.pcolormesh(mat,cmap=cm,norm=norm,edgecolors='white',linewidth=.55)
a.set_ylim(13,0);a.set_xticks(np.arange(6)+.5,['0','1','2','3','4','Ens.']);a.set_yticks(np.arange(13)+.5,[f'L{i:02d}' for i in range(1,13)]+['Macro']);a.tick_params(length=0,pad=3);a.set_xlabel('Checkpoint',labelpad=3)
a.axhline(12,color='#4B4F50',lw=.7);a.axvline(5,color='#4B4F50',lw=.7)
for i in range(13):
 for j in range(6):
  val=mat[i,j];lab=f'{val:.3f}' if val>=.001 else f'{val:.0e}'.replace('e-0','e-');rgb=cm(norm(val))[:3];lum=np.dot(rgb,[.299,.587,.114]);a.text(j+.5,i+.5,lab,ha='center',va='center',fontsize=7.5,color=INK if lum>.59 else 'white')
cbar=f.colorbar(im,cax=f.add_axes([.105,.517,.395,.009]),orientation='horizontal',ticks=[1e-5,1e-3,1e-1,1]);cbar.solids.set_rasterized(False);cbar.ax.minorticks_off();cbar.ax.set_xticklabels(['$10^{-5}$','$10^{-3}$','$10^{-1}$','1']);cbar.ax.tick_params(length=2,pad=2);cbar.outline.set_linewidth(.45);cbar.set_label('SD(predicted) / SD(observed)',fontsize=8,labelpad=2)
# B: raincloud distribution; actual data, macro summaries and zero counts retained.
b=f.add_axes([.66,.606,.315,.331]);rng=np.random.default_rng(5)
for j,key in enumerate(['atom','residue','image']):
 vals=np.array(d['branches'][key]['values']);nz=vals[vals>0];lv=np.log10(nz);yy=np.linspace(lv.min()-.15,lv.max()+.15,180);den=gaussian_kde(lv)(yy);ww=.29*den/den.max()
 b.fill_betweenx(10**yy,j-.08-ww,j-.08,facecolor=G,alpha=.20,edgecolor=G,lw=.7)
 b.scatter(j+rng.uniform(.02,.14,len(nz)),nz,s=7,facecolor='white',edgecolor=G,alpha=.68,lw=.5,zorder=2)
 b.scatter(np.full(5,j+.30),d['checkpoint_branch_macro'][j],s=11,marker='s',color=P,zorder=3)
 b.scatter(j-.015,d['branches'][key]['mean'],marker='D',s=28,color=DARK,edgecolor='white',lw=.4,zorder=4)
b.set_yscale('log');b.set_ylim(7e-9,.12);b.set_yticks([1e-8,1e-6,1e-4,1e-2]);b.set_xlim(-.5,2.5);b.set_xticks(range(3),['Atom','Residue','Image']);b.set_ylabel('|Shapley contribution| (kcal/mol)',labelpad=2);clean(b);b.yaxis.set_minor_locator(NullLocator());b.grid(axis='y',color='#ECEDED',lw=.5);b.set_axisbelow(True)
handles=[Line2D([],[],ls='',marker='o',ms=3.5,mfc='white',mec=G,mew=.6,label='Landscape'),Line2D([],[],ls='',marker='s',ms=3.5,color=P,label='Ckpt'),Line2D([],[],ls='',marker='D',ms=4,color=DARK,label='Ensemble')]
f.legend(handles=handles,loc='center',bbox_to_anchor=(.808,.565),ncol=3,frameon=False,fontsize=7.5,handlelength=.6,handletextpad=.35,columnspacing=.7)
total=sum(d['branches'][key]['mean'] for key in ['atom','residue','image']);txt(.575,.531,'Share',fontsize=7.5,color=GRAY);txt(.575,.508,'Zero',fontsize=7.5,color=GRAY)
for j,key in enumerate(['atom','residue','image']):
 xpos=.66+.315*(j+.5)/3
 txt(xpos,.531,f"{100*d['branches'][key]['mean']/total:.2f}%",ha='center',fontsize=8,color=DARK,weight='bold');txt(xpos,.508,f"{d['branches'][key]['zeros']}/46",ha='center',fontsize=7.5,color=GRAY)
# D: one landscape key shared by baseline, amplitude and error.
txt(.025,.444,'C',fontsize=11,weight='bold')
leg=[Patch(facecolor=CB1,edgecolor='#004A75',lw=.5,label='Mask only'),Patch(facecolor=CB2,edgecolor='#8A5A00',lw=.5,label='Complete alignment')]
f.legend(handles=leg,loc='upper right',bbox_to_anchor=(.985,.455),ncol=2,frameon=False,fontsize=7.5,handlelength=1.3,handletextpad=.5,columnspacing=1.3)
rect=[.338,.129,.29,.258];rect2=[.75,.129,.225,.258];left=f.add_axes(rect);right=f.add_axes(rect2,sharey=left);yy=np.arange(13)
for ax in [left,right]:
 ax.set_ylim(12.65,-.65);ax.axvline(0,color='#6B6D6E',lw=.65,zorder=1);ax.axhspan(11.6,12.65,color='#EFF2F4',zorder=0);ax.set_yticks([]);ax.spines[['left','right','top']].set_visible(False);ax.tick_params(length=2.3);ax.set_axisbelow(True);ax.grid(axis='x',color='#EEEFEF',lw=.5)
for j,c in enumerate(dd['conditions'][1:]):
 z=v[c];sd=np.r_[z['SD_change_percent'],z['macro_SD_change_percent']];ce=1000*np.r_[z['cMAE_change'],z['macro_cMAE_change']];off=[-.17,.17][j];style=dict(facecolor=CB1 if j==0 else CB2,edgecolor='#004A75' if j==0 else '#8A5A00',hatch=None,linewidth=.5)
 left.barh(yy+off,sd,height=.28,**style);right.barh(yy+off,ce,height=.28,**style)
left.set_xlim(-6.5,6.5);left.set_xticks([-6,-3,0,3,6]);left.set_xlabel('SD ratio change (%)',labelpad=3)
right.set_xlim(-.23,2.15);right.set_xticks([0,1,2]);right.set_xlabel('ΔcMAE (10⁻³ pKd)',labelpad=3)
txt(.218,.405,'Baseline SD ratio',ha='center',fontsize=7.5,color=GRAY)
orig=np.r_[v['original']['SD_ratio'],v['original']['macro_SD_ratio']]
for i,(name,val) in enumerate(zip(dd['landscapes']+['Macro'],orig)):
 y=rect[1]+rect[3]*(1-(i+.65)/13.3);weight='bold' if i==12 else 'normal';txt(.065,y,name,va='center',fontsize=7.5,weight=weight);lab=f'{val:.4f}' if val>.001 else f'{val:.1e}'.replace('e-0','e-');txt(.218,y,lab,ha='center',va='center',fontsize=7.5,color=DARK if i==12 else GRAY,weight=weight)
# Compact controls occupy a single aligned footer.
tab=f.add_axes([.065,.011,.91,.060]);tab.set_xlim(0,6);tab.set_ylim(0,3);tab.axis('off');tab.add_patch(Rectangle((0,0),6,2,facecolor='#F4F5F5',edgecolor='none'))
for x,label in [(3.35,'Original'),(4.45,'Mask'),(5.50,'Complete')]:tab.text(x,2.45,label,ha='center',va='center',fontsize=7.5,color=GRAY)
for y,label,key,den in [(1.5,'Equal residue vectors','equal_residue_vectors',1055),(.5,'Zero output response','output_zeros',211)]:
 tab.text(.05,y,f'{label} (/{den})',va='center',fontsize=7.5)
 for x,c in zip([3.35,4.45,5.50],dd['conditions']):tab.text(x,y,str(v[c][key]),ha='center',va='center',fontsize=8.5,weight='bold',color=DARK if v[c][key]==0 else INK)
f.savefig(O/'figure5.png',dpi=400)
sys.path.insert(0,str(Path(__file__).resolve().parent))
from visual_qa import audit_layout
issues=audit_layout(f);(O/'layout_audit.json').write_text(json.dumps(issues,indent=2,default=str));print(json.dumps(issues,default=str))
assert not issues,issues
f.savefig(O/'figure5.pdf');f.savefig(O/'figure5.svg');Image.open(O/'figure5.png').convert('L').save(O/'figure5_grayscale.png')
