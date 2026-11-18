"""Family-macro affinity and paired mutation evaluation; inputs are already qualified."""
import argparse, csv, json, math
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np

def finite(x):
 try: return math.isfinite(float(x))
 except (ValueError, TypeError): return False

def rank_correlation(a,b):
 def ranks(x):
  return np.array([sum(y<v for y in x)+(sum(y==v for y in x)+1)/2 for v in x],float)
 a,b=ranks(a),ranks(b)
 if len(a)<2 or np.ptp(a)==0 or np.ptp(b)==0: return None
 return float(np.corrcoef(a,b)[0,1])

def landscape_metrics(y,x,pos):
 y=np.asarray(y,float); x=np.asarray(x,float); pos=np.asarray(pos); counts=Counter(pos); w=np.array([1/len(counts)/counts[p] for p in pos]); q=x-y; qc=q-w@q; xc=x-w@x; yc=y-w@y
 vx=float(w@(xc*xc)); vy=float(w@(yc*yc)); cov=float(w@(xc*yc))
 return dict(ordinary_mae=float(w@abs(q)),cmae=float(w@abs(qc)),total_mse=float(w@(q*q)),mean_mse=float((w@q)**2),centered_mse=float(w@(qc*qc)),predicted_variance=vx,observed_variance=vy,covariance=cov,sd_ratio=math.sqrt(vx/vy) if vy>0 else None,amplitude_distance_from_one=abs(math.sqrt(vx/vy)-1) if vy>0 else None,position_mean_spearman=rank_correlation([float(np.mean(x[pos==p])) for p in sorted(counts)],[float(np.mean(y[pos==p])) for p in sorted(counts)]),zero_response=float(w@(x==0)))

def read(path):
 with open(path,encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
def mean(xs):
 x=[float(v) for v in xs if v is not None]; return float(np.mean(x)) if x else None

def run(config_path):
 config_path=Path(config_path).resolve(); c=json.loads(config_path.read_text(encoding='utf8')); base=config_path.parent
 def path(s): return (base / s).resolve()
 units=read(path(c['units'])); pairs=read(path(c['pairs'])); pred=read(path(c['predictions'])); U={u['unit_id']:u for u in units}; P={p['pair_id']:p for p in pairs}
 if len(U)!=len(units) or len(P)!=len(pairs): raise ValueError('duplicate unit or pair ID')
 methods=c['methods']; data={}; seeds={}; identities={}
 if not methods: raise ValueError('at least one method is required')
 for m,mc in methods.items():
  seeds[m]=[str(s) for s in mc['seeds']]; data[m]={}; identities[m]={}
  multiplier=mc.get('multiplier',1)
  if isinstance(multiplier,bool) or not isinstance(multiplier,(int,float)) or not finite(multiplier): raise ValueError(f'{m}: multiplier must be a finite JSON number')
  if mc['type'] not in ('endpoint','direct','score','zero') or not seeds[m] or not mc.get('identity'): raise ValueError('method requires valid type, seeds and identity')
  for row in pred:
   if row['method']!=m: continue
   k=(row['seed'],row['member_id'])
   if k in data[m]: raise ValueError('duplicate method, seed, member_id')
   if row['status']=='success' and finite(row['value']):
    converted=float(row['value'])*multiplier
    data[m][k]=converted if finite(converted) and abs(converted)<=math.sqrt(np.finfo(float).max)/4 else None; identities[m][k]=row.get('checkpoint','')
   else: data[m][k]=None
 def val(m,s,k): return data[m].get((s,k))
 def effect(m,s,p):
  mc=methods[m]; typ=mc['type']
  if typ=='zero': return 0.
  if typ in ('direct','score'): return val(m,s,p['pair_id'])
  a=val(m,s,p['wt_unit_id']); b=val(m,s,p['mutant_unit_id'])
  if a is None or b is None: return None
  ia=identities[m].get((s,p['wt_unit_id'])); ib=identities[m].get((s,p['mutant_unit_id']))
  if not ia or ia!=ib: raise ValueError('paired endpoints require the same recorded checkpoint or fixed ensemble identity')
  delta=b-a
  return delta if finite(delta) and abs(delta)<=math.sqrt(np.finfo(float).max)/4 else None
 def valid(m,p): return all(effect(m,s,p) is not None for s in seeds[m])
 def unitvalid(m,u): return methods[m]['type']=='endpoint' and all(val(m,s,u['unit_id']) is not None for s in seeds[m])
 for p in pairs:
  if p['wt_unit_id'] not in U or p['mutant_unit_id'] not in U: raise ValueError('pair has unknown endpoints')
  if p['family']!=U[p['wt_unit_id']]['family'] or p['family']!=U[p['mutant_unit_id']]['family']: raise ValueError('pair family disagrees with endpoints')
  if not finite(U[p['wt_unit_id']].get('observed')) or not finite(U[p['mutant_unit_id']].get('observed')): raise ValueError('paired endpoint labels must be finite')
  if not p['landscape'] or not p['position'] or not finite(p['observed_delta']): raise ValueError('pair requires landscape, position and finite observed_delta')
  if not math.isclose(float(p['observed_delta']),float(U[p['mutant_unit_id']]['observed'])-float(U[p['wt_unit_id']]['observed']),rel_tol=0,abs_tol=float(c.get('delta_tolerance',1e-8))): raise ValueError('observed_delta disagrees with endpoints')
 common=[p for p in pairs if all(valid(m,p) for m in methods)]
 lands=defaultdict(list)
 for p in common: lands[(p['family'],p['landscape'])].append(p)
 local=[p for pp in lands.values() if len(set(p['position'] for p in pp))>=int(c.get('min_positions',5)) for p in pp]
 result={'input':{'unit':c['unit'],'direction':c['direction'],'model_configurations':methods},'coverage':{'candidates':len(pairs),'successful_by_method':{m:sum(valid(m,p) for p in pairs) for m in methods},'common_pairs':len(common),'local_pairs':len(local),'local_landscapes':len({(p['family'],p['landscape']) for p in local}),'local_families':len({p['family'] for p in local})},'scopes':{}}
 for scope,ps in [('ordinary',common),('local',local)]:
  grouped=defaultdict(list)
  for p in ps: grouped[(p['family'],p['landscape'])].append(p)
  families=sorted({f for f,l in grouped}); out={'pairs':len(ps),'landscapes':len(grouped),'families':len(families),'methods':{},'contrasts':{}}; famvals={}
  for m,mc in methods.items():
   seed_fam=defaultdict(lambda:defaultdict(list)); defined=0; possible=0; direction=Counter()
   for s in seeds[m]:
    for (f,l),pp in sorted(grouped.items()):
     y=np.array([float(p['observed_delta']) for p in pp]); x=np.array([effect(m,s,p) for p in pp]); v=landscape_metrics(y,x,[p['position'] for p in pp]); possible+=1; defined+=v['position_mean_spearman'] is not None
     if mc['type']=='score':
      for k in ['ordinary_mae','cmae','total_mse','mean_mse','centered_mse','predicted_variance','observed_variance','covariance','sd_ratio','amplitude_distance_from_one']: v.pop(k)
     if mc['type']=='endpoint':
      w=np.array([1/len(set(p['position'] for p in pp))/sum(q['position']==p['position'] for q in pp) for p in pp]); ew=np.array([val(m,s,p['wt_unit_id'])-float(U[p['wt_unit_id']]['observed']) for p in pp]); em=np.array([val(m,s,p['mutant_unit_id'])-float(U[p['mutant_unit_id']]['observed']) for p in pp]); v.update(endpoint_mse=float(w@((ew*ew+em*em)/2)),shared_mse=float(w@(((ew+em)/2)**2)),difference_mse=float(w@((em-ew)**2)/4))
     for a,b in zip(y,x):
      direction['all']+=1; direction['three_class_correct']+=int(np.sign(a)==np.sign(b)); direction['experimental_nonzero']+=int(a!=0); direction['strict_correct']+=int(a!=0 and b!=0 and np.sign(a)==np.sign(b)); direction['both_nonzero']+=int(a!=0 and b!=0); direction['prediction_zero']+=int(b==0)
     for k,vv in v.items(): seed_fam[(f,s)][k].append(vv)
   family={f:{k:mean([mean(seed_fam[(f,s)].get(k,[])) for s in seeds[m]]) for k in {k for s in seeds[m] for k in seed_fam[(f,s)]}} for f in families}; famvals[m]=family
   keys={k for v in family.values() for k in v}; out['methods'][m]={'metrics':{k:mean([family[f].get(k) for f in families]) for k in sorted(keys)},'family_values':family,'rank_defined_landscape_seed':defined,'rank_possible_landscape_seed':possible,'direction_counts_across_seed_pair_incidences':dict(direction)}
  if families:
   rng=np.random.default_rng(int(c.get('bootstrap_seed_local' if scope=='local' else 'bootstrap_seed_ordinary',20260921 if scope=='local' else 20261021))); n=int(c.get('bootstrap_draws_local' if scope=='local' else 'bootstrap_draws_ordinary',10000 if scope=='local' else 2000)); draw=rng.integers(0,len(families),(n,len(families)))
   ref=c.get('reference')
   if ref in methods:
    for m in methods:
     if m==ref: continue
     contrast={}; commonkeys=set(out['methods'][m]['metrics']) & set(out['methods'][ref]['metrics'])
     for k in commonkeys:
      if k=='position_mean_spearman': continue
      a=[famvals[ref][f].get(k) for f in families]; b=[famvals[m][f].get(k) for f in families]
      if any(x is None for x in a+b): continue
      d=np.array(a)-np.array(b); contrast[k]={'reference_minus_method':float(d.mean()),'ci95':np.quantile(d[draw].mean(axis=1),[.025,.975]).tolist()}
     out['contrasts'][m]=contrast
   out['bootstrap']={'draws':n,'seed':int(c.get('bootstrap_seed_local' if scope=='local' else 'bootstrap_seed_ordinary',20260921 if scope=='local' else 20261021)),'family_order':families}
  # Restrict ordering comparisons to landscapes defined for every method and seed.
  rank_methods=c.get('rank_methods',[m for m in methods if methods[m]['type']!='zero']); rank_common=[]
  for (f,l),pp in sorted(grouped.items()):
   y=np.array([float(p['observed_delta']) for p in pp]); pos=[p['position'] for p in pp]
   if all(landscape_metrics(y,[effect(m,s,p) for p in pp],pos)['position_mean_spearman'] is not None for m in rank_methods for s in seeds[m]): rank_common.append((f,l))
  out['rank_common_landscapes']=len(rank_common); out['rank_common_families']=len({f for f,l in rank_common}); out['rank_common_metrics']={}
  rank_family={}
  for m in rank_methods:
   ff=defaultdict(list)
   for f,l in rank_common:
    pp=grouped[(f,l)]; y=np.array([float(p['observed_delta']) for p in pp]); pos=[p['position'] for p in pp]
    ff[f].append(mean([landscape_metrics(y,[effect(m,s,p) for p in pp],pos)['position_mean_spearman'] for s in seeds[m]]))
   rank_family[m]={f:mean(v) for f,v in ff.items()}
   out['rank_common_metrics'][m]=mean(list(rank_family[m].values()))
  rank_fams=sorted({f for f,l in rank_common}); ref=c.get('reference')
  if rank_fams and ref in rank_methods:
   seed=int(c.get('bootstrap_seed_local' if scope=='local' else 'bootstrap_seed_ordinary',20260921 if scope=='local' else 20261021)); n=int(c.get('bootstrap_draws_local' if scope=='local' else 'bootstrap_draws_ordinary',10000 if scope=='local' else 2000)); rank_draw=np.random.default_rng(seed).integers(0,len(rank_fams),(n,len(rank_fams)))
   for m in rank_methods:
    if m==ref: continue
    d=np.array([rank_family[ref][f]-rank_family[m][f] for f in rank_fams]); out['contrasts'][m]['position_mean_spearman']={'reference_minus_method':float(d.mean()),'ci95':np.quantile(d[rank_draw].mean(axis=1),[.025,.975]).tolist(),'landscapes':len(rank_common),'families':len(rank_fams),'scope':'all rank_methods common-defined','family_order':rank_fams}
  result['scopes'][scope]=out
 # Absolute common members use only endpoint methods; no invented absolute for direct or score.
 endpoint=[m for m in methods if methods[m]['type']=='endpoint']; eligible=[u for u in units if finite(u.get('observed')) and endpoint and all(unitvalid(m,u) for m in endpoint)]
 absout={'candidates':len(units),'successful_by_method':{m:sum(finite(u.get('observed')) and unitvalid(m,u) for u in units) for m in endpoint},'common_units':len(eligible),'families':len(set(u['family'] for u in eligible)),'methods':{},'contrasts':{}}
 for m in endpoint:
  family={f:mean([mean([abs(val(m,s,u['unit_id'])-float(u['observed'])) for u in eligible if u['family']==f]) for s in seeds[m]]) for f in sorted(set(u['family'] for u in eligible))}; absout['methods'][m]={'mae':mean(list(family.values())),'family_values':family}
 af=sorted({u['family'] for u in eligible})
 if af:
  n=int(c.get('bootstrap_draws_ordinary',2000)); seed=int(c.get('bootstrap_seed_ordinary',20261021)); draw=np.random.default_rng(seed).integers(0,len(af),(n,len(af)))
  for m in endpoint:
   a=np.array([absout['methods'][m]['family_values'][f] for f in af]); absout['methods'][m]['ci95']=np.quantile(a[draw].mean(axis=1),[.025,.975]).tolist()
  ref=c.get('reference')
  if ref in endpoint:
   for m in endpoint:
    if m!=ref:
     d=np.array([absout['methods'][ref]['family_values'][f]-absout['methods'][m]['family_values'][f] for f in af]); absout['contrasts'][m]={'reference_minus_method':float(d.mean()),'ci95':np.quantile(d[draw].mean(axis=1),[.025,.975]).tolist()}
  absout['bootstrap']={'draws':n,'seed':seed,'family_order':af}
 result['absolute']=absout
 # Linked absolute scopes deduplicate endpoint unit IDs independently for each scope.
 result['linked_absolute']={}
 for scope,ps in [('ordinary',common),('local',local)]:
  ids={p[k] for p in ps for k in ('wt_unit_id','mutant_unit_id')}; subset=[u for u in eligible if u['unit_id'] in ids]; bymethod={}
  for m in endpoint:
   fv={f:mean([mean([abs(val(m,s,u['unit_id'])-float(u['observed'])) for u in subset if u['family']==f]) for s in seeds[m]]) for f in sorted({u['family'] for u in subset})}; bymethod[m]={'mae':mean(list(fv.values())),'family_values':fv}
  result['linked_absolute'][scope]={'units':len(subset),'families':len({u['family'] for u in subset}),'methods':bymethod}
 return result

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('config'); ap.add_argument('--output',required=True); args=ap.parse_args(); result=run(args.config); out=Path(args.output); out.mkdir(parents=True,exist_ok=True); (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf8'); lines=['# Evaluation report','',f"Candidates: {result['coverage']['candidates']}; common pairs: {result['coverage']['common_pairs']}; local pairs: {result['coverage']['local_pairs']}.",f"Absolute common units: {result['absolute']['common_units']} in {result['absolute']['families']} families.",f"Absolute method estimates and intervals: {json.dumps({m:{k:v for k,v in d.items() if k != 'family_values'} for m,d in result['absolute']['methods'].items()})}; contrasts: {json.dumps(result['absolute']['contrasts'])}.",'','Weights: positions equally within landscape, landscapes equally within family, families equally. Seeds scored separately.','']
 for scope,v in result['scopes'].items():
  lines += [f"## {scope}: {v['pairs']} pairs, {v['landscapes']} landscapes, {v['families']} families",'']
  for m,d in v['methods'].items(): lines += [f"- {m}: {json.dumps(d['metrics'],ensure_ascii=False)}; rank coverage {d['rank_defined_landscape_seed']}/{d['rank_possible_landscape_seed']}"]
  lines += [f"Common defined ranks: {v['rank_common_landscapes']} landscapes in {v['rank_common_families']} families; common rank means: {json.dumps(v['rank_common_metrics'])}.", f"Reference-minus-method contrasts and 95% family bootstrap intervals: {json.dumps(v['contrasts'])}.", '']
 lines += ['Machine-readable summary contains family-level distributions, direction counts, linked-unit results and bootstrap metadata.']
 (out/'report.md').write_text('\n'.join(lines),encoding='utf8')
if __name__=='__main__': main()
