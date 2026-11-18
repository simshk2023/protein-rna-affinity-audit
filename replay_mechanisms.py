import argparse,json,re
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--inputs',required=True);ap.add_argument('--methods',required=True);ap.add_argument('--figures',required=True);ap.add_argument('--output',required=True);a=ap.parse_args();I=Path(a.inputs);O=Path(a.output);O.mkdir(parents=True,exist_ok=True);M=Path(a.methods);D=json.loads((M/'saved_inputs.json').read_text());P,U=D['P'],D['D'];native=json.loads((M/'m2_predictions.json').read_text())['pair_native_deltas'];rows=[]
for c,ck in native.items():
 for k,xx in ck.items():
  lands=defaultdict(list)
  for p,v in xx.items():
   if U[P[p]['mutant_unit_id']]['PDB'].upper()!='1JBR':lands[(P[p]['family_id'],P[p]['landscape_id'])].append(p)
  for sc in ['main','bshape']:
   ff=defaultdict(list);n=0
   for (f,l),pp in lands.items():
    pos=[re.fullmatch(r'[A-Z](\d+)[A-Z]',U[P[p]['mutant_unit_id']]['mutation_tokens'])[1] for p in pp];cnt=Counter(pos)
    if sc=='bshape' and len(cnt)<5:continue
    n+=len(pp);w=np.array([1/len(cnt)/cnt[v] for v in pos]);x=np.array([-xx[p]/1.3635606572452885 for p in pp]);y=np.array([float(P[p]['delta_pkd']) for p in pp]);q=x-y;ff[f].append([w@abs(q),w@abs(q-w@q)])
   vals=np.mean([np.mean(v,axis=0) for v in ff.values()],axis=0);rows.append(dict(condition=c,checkpoint=k,scope=sc,pairs=n,ordinary=float(vals[0]),centered=float(vals[1])))
(O/'m2_condition_metrics.json').write_text(json.dumps(rows,indent=2))
training=json.loads((I/'training.json').read_text());validation=json.loads((I/'validation.json').read_text());diag=json.loads((I/'diagnostics.json').read_text());fig=json.loads((Path(a.figures)/'figure7/plot_data.json').read_text());checks=[];fields={'absolute':'absolute_family_macro_mae','ordinary':'ordinary_delta_family_macro_mae','centered':'b_shape_cmae_family_macro'}
for sc,src in [('training',training),('validation',validation)]:
 for metric in fig['metrics']:
  for arm in fig['arms']:
   if metric in fields:values=[src['metrics'][f's{s}_{arm}'][fields[metric]] for s in fig['seeds']]
   else:values=diag['train_sp_amp_b' if sc=='training' else 'val_sp_amp_b'][arm][('spearman' if metric=='rank' else 'amplitude')+'_per_seed']
   expected=fig['matrix'][sc][metric][arm];delta=max(abs(x-y) for x,y in zip(values,expected['seeds']));delta=max(delta,abs(float(np.mean(values))-expected['mean']),abs(float(np.std(values,ddof=1))-expected['sd']));checks.append(dict(scope=sc,metric=metric,arm=arm,max_difference=delta))
assert max(x['max_difference'] for x in checks)<1e-11
(O/'three_arm_check.json').write_text(json.dumps(dict(status='PASS',cells=len(checks),seed_values=3*len(checks),max_difference=max(x['max_difference'] for x in checks)),indent=2));print(json.dumps({'three_arm_cells':len(checks),'m2_condition_scope_rows':len(rows)}))

ref=json.loads((I/'m2_reference.json').read_text());ds=[abs(x[k]-ref['metrics'][x['condition']][x['checkpoint']][x['scope']][{'ordinary':'ordinary_delta_pkd_mae','centered':'centered_delta_pkd_mae'}[k]]) for x in rows for k in ['ordinary','centered']];assert max(ds)<1e-11;(O/'m2_mechanism_check.json').write_text(json.dumps(dict(status='PASS',metrics=len(ds),max_difference=max(ds)),indent=2))
