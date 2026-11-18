"""CPU analysis of saved PRA formula terms. Standard library, no feature extraction."""
import json, math, statistics as st, hashlib, platform, socket, datetime
from pathlib import Path
from collections import defaultdict, Counter
import argparse
ap=argparse.ArgumentParser();ap.add_argument('--inputs',required=True);ap.add_argument('--output',required=True);ap.add_argument('--formula',required=True);a=ap.parse_args();O=Path(a.output);O.mkdir(parents=True,exist_ok=True);Q=O/'restricted';Q.mkdir(exist_ok=True);import shutil;shutil.copyfile(Path(a.inputs)/'pra_inputs.json',Q/'inputs.json');K=1.3635606572452885
data=json.loads((Q/'inputs.json').read_text()); pairs=data['pairs'];ends=data['endpoints']
FORM=json.loads(Path(a.formula).read_text(encoding='utf-8'))
assert set(FORM)=={'ss+str','ds+enz'},'PRA formula needs both branches'
for _b,_s in FORM.items():assert set(_s)=={'intercept','terms'},(_b,'formula schema')
def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False),encoding='utf-8')
def distribution(v):
 v=sorted(v)
 if not v:return {'n':0,'mean':None,'median':None,'min':None,'max':None,'q25':None,'q75':None}
 def quant(q):
  a=(len(v)-1)*q;i=int(a);return v[i]+(v[min(i+1,len(v)-1)]-v[i])*(a-i)
 return {'n':len(v),'mean':st.mean(v),'median':st.median(v),'min':v[0],'max':v[-1],'q25':quant(.25),'q75':quant(.75)}
def classify(c,d):
 if all(x==0 for x in c):return 'all_terms_unchanged'
 if d==0:return 'contributions_cancel'
 return 'nonzero_response'
assert classify([0,0],0)=='all_terms_unchanged'
assert classify([1,-1],0)=='contributions_cancel'
assert classify([1,-.5],.5)=='nonzero_response'
assert classify([.001],.001)=='nonzero_response' and float(f'{1.001:.2f}')==float(f'{1:.2f}')
formula_errors=[];rounded_errors=[];source_hash_match=[]
for eid,e in ends.items():
 f=e['features'];b=e['prediction']['branch'];spec=FORM[b]
 vals=[beta*f[k] for k,beta in spec['terms'].items()]
 g=vals[0]
 for x in vals[1:]:g+=x
 g+=spec['intercept']
 formula_errors.append(abs(g-e['prediction']['G_unrounded']))
 rounded_errors.append(float(f'{g:.2f}')==e['prediction']['G_2dp']==e['identity']['delta_g'])
 source_hash_match.append(e['identity']['src_sha256'].lower()=='667157288de06e6ef57c98f65c4c0ea2ec8be95a17d1c2cdb566046e79ebd88f')
 e['terms']={k:beta*f[k] for k,beta in spec['terms'].items()}
assert max(formula_errors)<1e-10 and all(rounded_errors) and all(source_hash_match)
records=[];sumerrors=[]
for p in pairs:
 w=ends[p['wt_endpoint_id']];m=ends[p['mut_endpoint_id']];b=p['branch'];spec=FORM[b]
 terms={k:{'wt':w['features'][k],'mut':m['features'][k],'difference':m['features'][k]-w['features'][k],'coefficient':beta,'contribution_G':beta*(m['features'][k]-w['features'][k]),'contribution_pKd':-beta*(m['features'][k]-w['features'][k])/K} for k,beta in spec['terms'].items()}
 c=[x['contribution_G'] for x in terms.values()];d=m['prediction']['G_unrounded']-w['prediction']['G_unrounded'];s=sum(c);den=sum(abs(x) for x in c)
 sumerrors.append(abs(s-d));cat=classify(c,d)
 records.append({**p,'terms':terms,'delta_G':d,'delta_pKd':-d/K,'contribution_sum_G':s,'category':cat,'sum_exact_zero':s==0,'near_zero_G_1e12':abs(d)<1e-12,'all_terms_near_zero_1e12':all(abs(x)<1e-12 for x in c),'rounded_zero':m['prediction']['G_2dp']==w['prediction']['G_2dp'],'rounding_only_zero':d!=0 and m['prediction']['G_2dp']==w['prediction']['G_2dp'],'cancellation':1-abs(s)/den if den else None})
assert max(sumerrors)<1e-10
def weights(rs):
 fam=defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
 for r in rs:fam[r['family_id']][r['landscape_id']][r['position']].append(r['pair_id'])
 return {p:1/len(fam)/len(ls)/len(ps)/len(ids) for ls in fam.values() for ps in ls.values() for ids in ps.values() for p in ids}
def moments(xs,ws):
 mu=sum(x*w for x,w in zip(xs,ws));return mu,sum(w*(x-mu)**2 for x,w in zip(xs,ws))
def variance_decomposition(values,rs,ws):
 # Weighted law of total variance and covariance: preserve negative allocations.
 keys=list(values[0]);land=defaultdict(list)
 for i,r in enumerate(rs):land[r['landscape_id']].append(i)
 means={}; masses={}
 for l,idx in land.items():
  mass=sum(ws[i] for i in idx);masses[l]=mass;means[l]={k:sum(ws[i]*values[i][k] for i in idx)/mass for k in keys}
 globalmean={k:sum(ws[i]*v[k] for i,v in enumerate(values)) for k in keys}
 total=[sum(v.values()) for v in values];gt=sum(ws[i]*x for i,x in enumerate(total))
 result={}
 for k in keys:
  bv=sum(masses[l]*(means[l][k]-globalmean[k])**2 for l in land)
  wv=sum(ws[i]*(values[i][k]-means[r['landscape_id']][k])**2 for i,r in enumerate(rs))
  bc=sum(masses[l]*(means[l][k]-globalmean[k])*(sum(means[l].values())-gt) for l in land)
  wc=sum(ws[i]*(values[i][k]-means[r['landscape_id']][k])*(total[i]-sum(means[r['landscape_id']].values())) for i,r in enumerate(rs))
  result[k]={'between_variance':bv,'within_variance':wv,'between_sd':math.sqrt(bv),'within_sd':math.sqrt(wv),'between_output_covariance':bc,'within_output_covariance':wc}
 tv=sum(ws[i]*(x-gt)**2 for i,x in enumerate(total));bv=sum(masses[l]*(sum(means[l].values())-gt)**2 for l in land);wv=tv-bv
 err=max(abs(sum(x['between_output_covariance'] for x in result.values())-bv),abs(sum(x['within_output_covariance'] for x in result.values())-wv))
 assert err<1e-10
 return {'terms':result,'output_variance':tv,'output_between_variance':bv,'output_within_variance':wv,'covariance_identity_error':err,'n_observations':len(rs),'n_landscapes':len(land)}
def scope(rs):
 cats=Counter(r['category'] for r in rs);ww=weights(rs);assert abs(sum(ww.values())-1)<1e-12
 return {'pairs':len(rs),'positions':len({(r['landscape_id'],r['position']) for r in rs}),'landscapes':len({r['landscape_id'] for r in rs}),'families':len({r['family_id'] for r in rs}),'units':len({r[k] for r in rs for k in ['wt_unit_id','mutant_unit_id']}),'categories':dict(cats),'exact_zero':sum(r['delta_G']==0 for r in rs),'near_zero_G_1e12':sum(r['near_zero_G_1e12'] for r in rs),'rounded_zero':sum(r['rounded_zero'] for r in rs),'rounding_only_zero':sum(r['rounding_only_zero'] for r in rs),'cancellation':distribution([r['cancellation'] for r in rs if r['cancellation'] is not None]),'hierarchical_category_weight':{c:sum(ww[r['pair_id']] for r in rs if r['category']==c) for c in cats},'sum_weights':sum(ww.values())}
main=[r for r in records if r['scope']=='main'];diag=[r for r in records if r['scope']=='diagnostic']
summary={'main':scope(main),'diagnostic_1JBR':scope(diag),'branches':{},'formulas':FORM,'weighting':'Within each branch: equal families, equal landscapes within family, equal mutation positions within landscape, equal pairs within position. Pair counts and raw distributions separately reported. Absolute scope deduplicates unit IDs, then equal families/landscapes/units.','energy_unit':'kcal/mol','pKd_conversion':-1/K,'near_zero_threshold_unit':'1e-12 kcal/mol','execution':{'hostname':socket.gethostname(),'platform':platform.platform(),'python':platform.python_version(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256((Q/'inputs.json').read_bytes()).hexdigest()}}
for b,spec in FORM.items():
 rs=[r for r in main if r['branch']==b];ww=weights(rs);ws=[ww[r['pair_id']] for r in rs];stats={}
 for k in spec['terms']:
  cs=[r['terms'][k]['contribution_G'] for r in rs]
  stats[k]={'changed_pairs':sum(r['terms'][k]['difference']!=0 for r in rs),'pairs':len(rs),'change_rate':sum(r['terms'][k]['difference']!=0 for r in rs)/len(rs),'hierarchical_change_rate':sum(ws[i] for i,r in enumerate(rs) if r['terms'][k]['difference']!=0),'absolute_contribution_G':distribution([abs(c) for c in cs]),'signed_contribution_G':distribution(cs),'hierarchical_mean_abs_contribution_G':sum(w*abs(c) for w,c in zip(ws,cs)),'hierarchical_mean_signed_contribution_G':sum(w*c for w,c in zip(ws,cs))}
 # Both endpoints enter the matched decomposition with half each pair weight.
 vr=[];vv=[];vw=[]
 for r,w in zip(rs,ws):
  for side in ['wt','mut']:
   vr.append(r);vv.append({k:spec['terms'][k]*r['terms'][k][side] for k in spec['terms']});vw.append(w/2)
 matched=variance_decomposition(vv,vr,vw)
 # WT-deduplicated absolute unit scope: shared WT only once.
 unique={}
 for r in rs:
  for side,uidkey in [('wt','wt_unit_id'),('mut','mutant_unit_id')]:
   uid=r[uidkey];v={k:spec['terms'][k]*r['terms'][k][side] for k in spec['terms']}
   if uid in unique:assert unique[uid][1]==v
   unique[uid]=(r,v)
 fam=defaultdict(lambda:defaultdict(list))
 for uid,(r,v) in unique.items():fam[r['family_id']][r['landscape_id']].append(uid)
 uw={u:1/len(fam)/len(ls)/len(us) for ls in fam.values() for us in ls.values() for u in us}
 absolute=variance_decomposition([v for r,v in unique.values()],[r for r,v in unique.values()],[uw[u] for u in unique])
 summary['branches'][b]={'response':scope(rs),'terms':stats,'matched_endpoint_variance':matched,'deduplicated_absolute_variance':absolute,'RNA_geometry_joint_unchanged_pairs':sum(all(r['terms'][k]['difference']==0 for k in ['Incl1','Buckle1','rise1']) for r in rs) if b=='ds+enz' else None}
 dump(Q/('weights_'+b.replace('+','_')+'.json'),ww)
assert (len(main),len(diag))==(106,11)
assert summary['main']['exact_zero']==42 and summary['main']['rounded_zero']==55
assert summary['main']['units']==121 and summary['main']['landscapes']==15 and summary['main']['families']==14
assert [len([r for r in main if r['branch']==b]) for b in FORM]==[13,93]
checks={'scope_106_plus_11':True,'endpoint_formula_1e10':True,'identity_sha_chain_success':True,'rounded_endpoint_agreement':True,'contribution_sum_1e10':True,'anchors_42_55':True,'hierarchical_weights':True,'variance_covariance_identity':True,'zero_denominator_is_null':all(r['cancellation'] is None for r in records if r['category']=='all_terms_unchanged'),'synthetic_classification_cases':True}
dump(O/'analysis.json',summary);dump(Q/'pair_contributions.json',records);dump(Q/'endpoint_terms.json',ends)
dump(O/'verification.json',{'all_pass':all(checks.values()),'checks':checks,'endpoints':len(ends),'max_formula_error_G':max(formula_errors),'max_contribution_sum_error_G':max(sumerrors),'execution':summary['execution']})
print(json.dumps({'all_pass':True,'main':summary['main'],'execution':summary['execution']},indent=2))
