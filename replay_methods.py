"""Method-specific saved-input aggregation, adapted from Hfq reconciliation, no inference."""
from pathlib import Path
import argparse,json,csv,collections,re
import numpy as np
from scipy.stats import spearmanr
ap=argparse.ArgumentParser();ap.add_argument('--inputs',required=True);ap.add_argument('--output',required=True);ap.add_argument('--service-predictions',required=True);a=ap.parse_args();I=Path(a.inputs);O=Path(a.output);O.mkdir(parents=True,exist_ok=True)
x=json.loads((I/'saved_inputs.json').read_text());P,D,U=x['P'],x['D'],x['U'];K=1.3635606572452885
ex={k:set(v) for k,v in x['ex'].items()}
meta=Path(__file__).resolve().parent/'metadata';drop=json.loads((meta/'source_alias_exclusions.json').read_text());dp=set(drop['pair_ids']);du=set(drop['unit_ids'])
def position(p):return int(re.fullmatch(r'[A-Z](\d+)[A-Z]',D[P[p]['mutant_unit_id']]['mutation_tokens'])[1])
def groups(ids):
 g=collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(list)))
 for p in ids:g[P[p]['family_id']][P[p]['landscape_id']][position(p)].append(p)
 return g
def weights(ids):
 g=groups(ids);w={}
 for ls in g.values():
  for ps in ls.values():
   for pp in ps.values():
    for p in pp:w[p]=1/len(g)/len(ls)/len(ps)/len(pp)
 return np.array([w[p] for p in ids])
def shape(ids):return sorted(p for ls in groups(ids).values() for ps in ls.values() if len(ps)>=5 for pp in ps.values() for p in pp)
def covered(m):return {p for p,r in P.items() if all(r['wt_unit_id'] in u and r['mutant_unit_id'] in u for u in U[m].values())}
ns={'spearmanr':spearmanr}
def score(ids,preds,units=None):
 ids=sorted(ids);bs=shape(ids);out={'pairs':len(ids),'families':len(groups(ids)),'bshape_pairs':len(bs),'bshape_landscapes':sum(len(v) for v in groups(bs).values()),'methods':{}}
 for m,ss in preds.items():
  records=[]
  for s,xall in ss.items():
   def err(pp,cent=False):
    y=np.array([float(P[p]['delta_pkd']) for p in pp]);x=np.array([xall[p] for p in pp]);q=x-y
    if cent:
     for ls in groups(pp).values():
      for ps in ls.values():
       member={p for v in ps.values() for p in v};pp2=[p for p in pp if p in member];ix=np.array([p in member for p in pp]);q[ix]-=weights(pp2)@q[ix]
    return float(weights(pp)@abs(q))
   r={'ordinary':err(ids),'local_ordinary':err(bs),'cmae':err(bs,True)};sr=[];rr=[]
   for ls in groups(bs).values():
    for ps in ls.values():
     pp=[p for v in ps.values() for p in v];w=weights(pp);y=np.array([float(P[p]['delta_pkd']) for p in pp]);x=np.array([xall[p] for p in pp]);sr.append(float(np.sqrt(w@(x-w@x)**2)/np.sqrt(w@(y-w@y)**2)))
     a=[np.mean([xall[p] for p in v]) for v in ps.values()];b=[np.mean([float(P[p]['delta_pkd']) for p in v]) for v in ps.values()]
     if len(set(a))>1:rr.append(float(ns['spearmanr'](a,b).statistic))
   r.update(sd_ratio=float(np.mean(sr)),sd_median=float(np.median(sr)),rank=float(np.mean(rr)) if rr else None)
   y=np.array([float(P[p]['delta_pkd']) for p in ids]);x=np.array([xall[p] for p in ids]);nz=y!=0;pos=y>0;neg=y<0
   r.update(direction=float(np.mean((np.sign(x[nz])==np.sign(y[nz])))),weak_fraction=float(np.mean(y[nz]<0)),sensitivity=float(np.mean(x[pos]>0)),specificity=float(np.mean(x[neg]<0)),strengthening_hits=int(sum(pos&(x>0))),weakening_hits=int(sum(neg&(x<0))),direction_n=int(sum(nz)),observed_strong=int(sum(pos)),observed_weak=int(sum(neg)))
   if units is not None and m!='zero':
    f=collections.defaultdict(list)
    for u in units:f[D[u]['family_id']].append(abs(U[m][s][u]-float(D[u]['label'])))
    r['absolute']=float(np.mean([np.mean(v) for v in f.values()]))
   records.append(r)
  out['methods'][m]={k:float(np.mean([r[k] for r in records])) if records[0][k] is not None else None for k in records[0]}
 if units is not None:out['units']=len(units)
 return out
def endpoint_scope(ids,methods,units=None):
 if units is None:units={P[p][s] for p in ids for s in ['wt_unit_id','mutant_unit_id']}
 pred={m:{s:{p:v[P[p]['mutant_unit_id']]-v[P[p]['wt_unit_id']] for p in ids} for s,v in U[m].items()} for m in methods};pred['zero']={'fixed':{p:0 for p in ids}}
 return score(ids,pred,units)

base=covered('copra')&covered('baseline')-ex['CENSORED_PAIRS']-dp
controllers={p for p in base if D[P[p]['mutant_unit_id']]['PDB'].upper()=='1JBR'}
def units_for(ms):return set.intersection(*(set(v) for m in ms for v in U[m].values()))-ex['CENSORED_UNITS']-du
scopes={}
for name,ms in [('S3',['copra','baseline','deeppnap','seq_baseline','omni_official','omni_controlled']),('S4',['copra','baseline','omni_official','omni_controlled']),('S6',['copra','baseline','prdelta']),('S8',['copra','baseline','svsbi','gbnl','prorb'])]:
 ids=base&set.intersection(*(covered(m) for m in ms))
 if name=='S6':ids-=controllers;units={P[p][k] for p in ids for k in ['wt_unit_id','mutant_unit_id']}
 else:units=units_for(ms)
 scopes[name]=endpoint_scope(ids,ms,units)
pra=json.loads((I/'pra_inputs.json').read_text());ids=base&{p['pair_id'] for p in pra['pairs']}-controllers;scopes['S9']=endpoint_scope(ids,['pra','copra','baseline'])
m2=json.loads((I/'m2_predictions.json').read_text())['endpoints']['all_indices_aligned']['ensemble'];man=json.loads((I/'m2_cohort_manifest.json').read_text());up={};ids=set()
for p in man['pairs']:
 if p['wt_endpoint_id'] in m2 and p['mut_endpoint_id'] in m2:
  ids.add(p['pair_id'])
  for side,ek in [('wt_unit_id','wt_endpoint_id'),('mutant_unit_id','mut_endpoint_id')]:up[P[p['pair_id']][side]]=-m2[p[ek]]/K
U['m2_aligned']={'ensemble':up};scopes['S10']=endpoint_scope(ids&base-controllers,['m2_aligned','copra','baseline'])
service=json.loads(Path(a.service_predictions).read_text());M=service['mcsm'];ids=set(M)-controllers-dp
pred={'mcsm':{'fixed':M},**{m:{s:{p:v[P[p]['mutant_unit_id']]-v[P[p]['wt_unit_id']] for p in ids} for s,v in U[m].items()} for m in ['copra','baseline']},'zero':{'fixed':{p:0 for p in ids}}};scopes['S7']=score(ids,pred)
pan=service['panther'];pids=sorted(pan);lands={lid:sorted(p for p in pids if P[p]['landscape_id']==lid) for lid in {P[p]['landscape_id'] for p in pids}};eligible=[pp for pp in lands.values() if len({position(p) for p in pp})>=5]
px={'panther':{'fixed':pan},**{m:{seed:{p:v[P[p]['mutant_unit_id']]-v[P[p]['wt_unit_id']] for p in pids} for seed,v in U[m].items()} for m in ['copra','baseline']}}
def rank(pp,x):
 ps=sorted({position(p) for p in pp});xx=[np.mean([x[p] for p in pp if position(p)==v]) for v in ps];yy=[np.mean([float(P[p]['delta_pkd']) for p in pp if position(p)==v]) for v in ps]
 return float(spearmanr(xx,yy).statistic) if len(set(xx))>1 and len(set(yy))>1 else None
common=[pp for pp in eligible if all(rank(pp,x) is not None for mm in px.values() for x in mm.values())];assert sum(map(len,common))==31
panout={}
for m,ss in px.items():
 seeds=[]
 for x in ss.values():
  fam=collections.defaultdict(list)
  for pp in common:fam[P[pp[0]]['family_id']].append(rank(pp,x))
  seeds.append(dict(zeros=sum(v==0 for v in x.values()),strict_direction=100*sum(np.sign(x[p])==np.sign(float(P[p]['delta_pkd'])) and x[p]!=0 and float(P[p]['delta_pkd'])!=0 for p in pids)/len(pids),rank=float(np.mean([np.mean(v) for v in fam.values()]))))
 panout[m]={k:float(np.mean([v[k] for v in seeds])) for k in seeds[0]}
scopes['S11']={'pairs':87,'rank_pairs':31,'rank_landscapes':len(common),'methods':panout}
(O/'method_results.json').write_text(json.dumps(scopes,indent=2,allow_nan=False));print(json.dumps({k:{f:v.get(f) for f in ['units','pairs','bshape_pairs']} for k,v in scopes.items()},indent=2))
