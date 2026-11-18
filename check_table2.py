import argparse,json,re,csv,sys
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('summary');ap.add_argument('--manuscript',required=True);ap.add_argument('--output',required=True);a=ap.parse_args();r=json.loads(Path(a.summary).read_text());tex=Path(a.manuscript).read_text(encoding='utf-8-sig');o=Path(a.output);o.mkdir(parents=True,exist_ok=True)
b=tex[tex.index('Measure & CoPRA'):tex.index('Sources: Tables S17')];lines=[x.strip() for x in b.splitlines() if '&' in x and not x.startswith('Measure')];rows=[]
def metrics(scope,k):return [r['scopes'][scope]['methods'][m]['metrics'][k] for m in ['copra','baseline','zero']]
def contrast(scope,k,invert=False):
 x=r['scopes'][scope]['contrasts']['copra'][k];v=[x['reference_minus_method'],*x['ci95']];return [-v[0],-v[2],-v[1]] if invert else v
specs={
'Absolute MAE, 752 units':(0,'absolute','mae'),
'Linked absolute MAE, 338-pair endpoints':(1,'ordinary','linked_absolute_mae'),
'Ordinary $\\Delta$MAE, 338 pairs':(2,'ordinary','ordinary_mae'),
'Linked absolute MAE, 241 units':(3,'local','linked_absolute_mae'),
'Ordinary $\\Delta$MAE':(4,'local','ordinary_mae'),
'cMAE':(5,'local','cmae'),
'Position-mean Spearman':(6,'local','position_mean_spearman'),
'Predicted / observed SD':(7,'local','sd_ratio'),
'$|\\mathrm{SD\\ ratio}-1|$':(8,'local','amplitude_distance_from_one')}
assert len(lines)==len(specs)
for line in lines:
 name=line.split('&')[0].strip();i,scope,metric=specs[name];vals=[]
 if i==0:vals=[r['absolute']['methods'][m]['mae'] for m in ['copra','baseline']]
 elif i in [1,3]:
  sc='ordinary' if i==1 else 'local';z=r['linked_absolute'][sc];vals=[z['methods'][m]['mae'] for m in ['copra','baseline']]
  if i==3:
   fs=sorted(z['methods']['copra']['family_values']);d=np.array([z['methods']['baseline']['family_values'][f]-z['methods']['copra']['family_values'][f] for f in fs]);dr=np.random.default_rng(20260921).integers(0,len(fs),(10000,len(fs)));vals += [float(d.mean()),*np.quantile(d[dr].mean(axis=1),[.025,.975]).tolist()]
 else:
  sc='ordinary' if i==2 else 'local';k={2:'ordinary_mae',4:'ordinary_mae',5:'cmae',6:'position_mean_spearman',7:'sd_ratio',8:'amplitude_distance_from_one'}[i];vals=metrics(sc,k)
  if k=='position_mean_spearman':vals=vals[:2]
  vals+=contrast(sc,k,k in ['position_mean_spearman','sd_ratio'])
 cells=line.split('&')[1:];expected=[]
 for c in cells:
  c=c.replace('$-$','-');expected += [float(x) for x in re.findall(r'(?<![A-Za-z])[-+]?\d+(?:\.\d+)?',c)]
 assert len(expected)==len(vals),(name,expected,vals)
 diffs=[abs(float(format(v,'.4f'))-e) for v,e in zip(vals,expected)];rows.append(dict(key=f'{scope}/{metric}/copra,baseline,zero,paired_effect',expected=expected,actual=vals,rounded_max_difference=max(diffs),status='PASS' if max(diffs)<1e-10 else 'FAIL'))
(o/'table2_keyed_check.json').write_text(json.dumps(rows,indent=2));assert all(x['status']=='PASS' for x in rows),[(x['key'],x['expected'],x['actual']) for x in rows if x['status']!='PASS'];print(json.dumps({'status':'PASS','rows':len(rows),'numeric_cells':sum(len(x['actual']) for x in rows)}))
