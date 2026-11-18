"""Three-panel revision: performance, joint family gains, landscape response."""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize,LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.transforms import Bbox
from PIL import Image
O=Path(__file__).resolve().parent;D=json.loads((O/'plot_data.json').read_text())
sys.path.insert(0,str(Path(__file__).resolve().parent))
from visual_qa import audit_layout
COL={'A':'#A81024','B':'#4C9C93','C':'#006176'};MARK={'A':'o','B':'s','C':'^'};INK='#252A2D';GRAY='#687479'
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none','text.color':INK,'axes.labelcolor':INK})
fig=plt.figure(figsize=(178/25.4,135/25.4),facecolor='white')
def text(x,y,t,size=8,**kw):return fig.text(x,y,t,fontsize=size,**kw)
def title(x,y,l,t):text(x,y,l,11,weight='bold');text(x+.036,y,t,9,weight='bold')
def clean(a):a.spines[['top','right']].set_visible(False);a.tick_params(length=2.5,pad=3)
title(.025,.969,'A','Performance')
X,Y,W,H=.147,.568,.813,.317
ax=fig.add_axes([X,Y,W,H]);ax.set_xlim(0,5);ax.set_ylim(6,0);ax.axis('off')
cm=LinearSegmentedColormap.from_list('performance',['#A81024','#FFF4C9','#006176'])
norms={'error':Normalize(-60,60),'rank':Normalize(-1,1),'amplitude':Normalize(0,2)}
for i,(s,a) in enumerate([(s,a) for s in ['training','validation'] for a in D['arms']]):
 for j,m in enumerate(D['metrics']):
  q=D['matrix'][s][m][a];key=m if m in ['rank','amplitude'] else 'error';v=norms[key](q['color_value'])
  ax.add_patch(Rectangle((j,i),1,1,facecolor=cm(v),edgecolor='white',lw=2))
  ax.text(j+.5,i+.5,f"{q['mean']:.3f}",fontsize=10,ha='center',va='center',color='white' if abs(v-.5)>.31 else INK)
 ax.text(-.09,i+.5,a,color=COL[a],fontsize=9,weight='bold',ha='right',va='center')
ax.axhline(3,color='white',lw=5)
for yy,t in [(.806,'Training'),(.648,'Validation')]:text(.035,yy,t,8.5,rotation=90,ha='center',va='center')
for j,t in enumerate(['Absolute MAE ↓','Ordinary ΔMAE ↓','Local cMAE ↓','Position Spearman ↑','Local SD ratio']):text(X+(j+.5)*W/5,.903,t,8,ha='center')
for xx,ww,key,ticks,t in [(X,.44,'error',[-60,0,60],'Error reduction vs A (%)'),(.653,.137,'rank',[-1,0,1],'Spearman'),(.823,.137,'amplitude',[0,1,2],'SD ratio')]:
 cb=fig.colorbar(ScalarMappable(norm=norms[key],cmap=cm),cax=fig.add_axes([xx,.510,ww,.011]),orientation='horizontal',ticks=ticks)
 cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face');cb.outline.set_visible(False);cb.ax.tick_params(length=2,labelsize=7.5,pad=2)
 text(xx+ww/2,.535,t,7.5,ha='center')
# B: merge two metrics onto one common gain axis, preserving each seed.
title(.025,.441,'B','Family-specific gains')
F=pd.read_csv(O/'family_seed_metrics.csv')
G=F.pivot(index=['family','seed','metric'],columns='arm',values='value').reset_index()
G['gain']=G.B-G.C
G.to_csv(O/'merged_gain_data.csv',index=False)
METCOL={'ordinary':'#006176','centered':'#C47781'}
METMARK={'ordinary':'o','centered':'s'}
b=fig.add_axes([.105,.130,.290,.254]);clean(b)
b.spines['left'].set_visible(False)
b.set_xlim(-.035,.090);b.set_ylim(-.48,3.48)
b.set_xticks([-.03,0,.03,.06,.09],['−0.03','0','0.03','0.06','0.09'])
b.set_yticks([3,2,1,0],['F1','F2','F3','F4']);b.tick_params(axis='y',length=0,pad=8)
for label in b.get_yticklabels():label.set_weight('bold')
b.axvline(0,color='#9DA9AE',lw=.8,zorder=0)
for yy in [.5,1.5,2.5]:b.axhline(yy,color='#E9EDEE',lw=.55,zorder=0)
for i,fam in enumerate(['F1','F2','F3','F4']):
 for j,metric in enumerate(['ordinary','centered']):
  z=G[(G.family==fam)&(G.metric==metric)].sort_values('seed')
  values=z.gain.to_numpy();mean=values.mean();yy=3-i+(.20 if j==0 else -.20)
  color=METCOL[metric];marker=METMARK[metric]
  # A thin data-range strip shows the actual minimum/maximum of three seeds.
  b.add_patch(Rectangle((values.min(),yy+.005),values.max()-values.min(),.14,facecolor=color,alpha=.10,edgecolor='none',zorder=1))
  b.scatter(values,yy+np.array([.025,.075,.125]),s=12,marker=marker,facecolor=color,edgecolor='white',lw=.3,alpha=.68,zorder=3)
  b.plot([0,mean],[yy-.07,yy-.07],color=color,lw=2.2,alpha=.7,solid_capstyle='round',zorder=2)
  b.scatter(mean,yy-.07,s=38,marker=marker,facecolor=color,edgecolor='white',lw=.6,zorder=4)
  b.text(1.075,(yy-.07+.48)/3.96,f'{mean:+.4f}',transform=b.transAxes,color=color,fontsize=7.5,ha='left',va='center',weight='bold',clip_on=False)
b.set_xlabel('Error reduction, B − C (pKd)',labelpad=5)
text(.418,.402,'Mean',7.5,color=GRAY)
for x,metric,label in [(.109,'ordinary','Ordinary ΔMAE'),(.281,'centered','Local cMAE')]:
 fig.add_artist(Line2D([x],[.409],transform=fig.transFigure,marker=METMARK[metric],color=METCOL[metric],linestyle='none',markersize=3.5))
 text(x+.012,.402,label,7.5,color=METCOL[metric])
# C: individual landscapes reveal heterogeneity masked by pooled summaries.
title(.501,.441,'C','Local response across landscapes')
L=pd.DataFrame(D['landscape_means']);S=pd.read_csv(O/'landscape_seed_metrics.csv');family_map={r['family_id']:r['family'] for r in D['family_counts']}
axes=[]
for k,scope in enumerate(['training','validation']):
 a=fig.add_axes([.591+k*.209,.130,.170,.254]);axes.append(a);clean(a)
 a.set_xlim(-.05,2.05);a.set_ylim(-.5,1.05);a.set_xticks([0,1,2]);a.set_yticks([-.5,0,.5,1]);a.axvline(1,color='#A4AFB2',lw=.7,ls='--');a.axhline(0,color='#A4AFB2',lw=.6)
 if k:a.set_yticklabels([])
 else:a.set_ylabel('Position Spearman',labelpad=3)
 a.set_title(('Training · 15 landscapes' if k==0 else 'Validation · 4 landscapes'),fontsize=7.5,pad=8)
 for lid,z in L[L.scope==scope].groupby('landscape_id'):
  z=z.set_index('arm').loc[D['arms']]
  # Connect matched landscape evaluations across objectives, never across splits.
  a.plot(z.amplitude,z['rank'],color='#A7B4B7',lw=.65,alpha=.7,zorder=1)
  if scope=='validation':
   fam=family_map[z.family_id.iloc[0]];mx=z.amplitude.mean();my=z['rank'].mean()
   labels={'F1':(.83,.32),'F2':(.22,-.29),'F3':(.64,.47),'F4':(.85,.08)}
   a.annotate(fam,(mx,my),xytext=labels[fam],textcoords='data',fontsize=7.5,color=GRAY,arrowprops=dict(arrowstyle='-',lw=.5,color=GRAY))
 for arm in D['arms']:
  ss=S[(S.scope==scope)&(S.arm==arm)]
  a.scatter(ss.amplitude,ss['rank'],marker=MARK[arm],s=5,color=COL[arm],alpha=.22,edgecolor='none',zorder=2)
  z=L[(L.scope==scope)&(L.arm==arm)]
  a.scatter(z.amplitude,z['rank'],marker=MARK[arm],s=20,facecolor=COL[arm],edgecolor='white',lw=.5,zorder=3)
  q=D['matrix'][scope]
  a.scatter(q['amplitude'][arm]['mean'],q['rank'][arm]['mean'],marker='D',s=40,facecolor='white',edgecolor=COL[arm],lw=1.1,zorder=5)
text(.779,.067,'Local SD ratio',8,ha='center')
# Single shared arm key. All methodological text is in the caption.
handles=[Line2D([],[],marker=MARK[a],color=COL[a],linestyle='none',markersize=4,label=t) for a,t in zip('ABC',['A  Absolute','B  + Endpoint','C  + Difference'])]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.004),ncol=3,frameon=False,fontsize=8,handletextpad=.5,columnspacing=2)
fig.canvas.draw();(O/'layout_audit.json').write_text(json.dumps(audit_layout(fig),indent=2,default=str))
for ext in ['png','pdf','svg']:fig.savefig(O/f'figure6.{ext}',dpi=400,facecolor='white')
Image.open(O/'figure6.png').convert('L').save(O/'figure6_grayscale.png',dpi=(400,400))
w,h=fig.get_size_inches()
fig.savefig(O/'panel_B.png',dpi=500,bbox_inches=Bbox.from_extents(.02*w,.038*h,.495*w,.475*h),facecolor='white')
print('Rendered revised three-panel figure, 178 × 143 mm.')
