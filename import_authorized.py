"""Reconstruct current inputs from provider CSV + curated nonnumeric identities + saved predictions."""
import argparse,csv,json,math,statistics,sys
from pathlib import Path
from collections import defaultdict
ap=argparse.ArgumentParser();ap.add_argument('--pronab',required=True);ap.add_argument('--predictions',required=True);ap.add_argument('--output',required=True);a=ap.parse_args()
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE/'toolkit'));import evaluate
out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
def rd(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def wr(n,rr):
 with (out/n).open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
meta=HERE/'metadata';raw=rd(a.pronab);units=rd(meta/'units.csv');pairs=rd(meta/'pairs.csv');pred=rd(a.predictions);ex=json.loads((meta/'censor_exclusions.json').read_text());aliases=json.loads((meta/'source_alias_exclusions.json').read_text());excluded=set(ex['CENSORED_UNITS'])|set(aliases['unit_ids'])
keep={r['unit_id'] for r in units}-excluded;links=defaultdict(list)
for r in rd(meta/'source_links.csv'):
 if r['unit_id'] not in keep:continue
 source=raw[int(r['source_row_number'])-1]
 if source['Entry_id']!=r['entry_id']:raise ValueError('source row/Entry_id mismatch; use recorded snapshot or explicitly re-align by Entry_id')
 text=source[r['raw_column']].strip()
 if any(x in text for x in '<>≤≥~'):raise ValueError('unexpected bound in retained member')
 v=float(text)
 if not math.isfinite(v) or v<=0:raise ValueError('invalid retained Kd')
 links[r['unit_id']].append(-math.log10(v))
u=[dict(unit_id=r['unit_id'],family=r['family_id'],observed=statistics.median(links[r['unit_id']])) for r in units if r['unit_id'] in keep];U={r['unit_id']:r for r in u};positions={r['pair_id']:r['position'] for r in rd(meta/'candidate_membership.csv')};eligible={r['pair_id'] for r in rd(meta/'candidate_membership.csv') if r['reason']!='censored_endpoint'}-set(aliases['pair_ids'])
p=[dict(pair_id=r['pair_id'],family=r['family_id'],landscape=r['landscape_id'],position=positions[r['pair_id']],wt_unit_id=r['wt_unit_id'],mutant_unit_id=r['mutant_unit_id'],observed_delta=U[r['mutant_unit_id']]['observed']-U[r['wt_unit_id']]['observed']) for r in pairs if r['pair_id'] in eligible]
wr('units.csv',u);wr('pairs.csv',p);wr('predictions.csv',pred)
c=dict(units='units.csv',pairs='pairs.csv',predictions='predictions.csv',unit='pKd',direction='mutant minus WT; higher is stronger',reference='baseline',methods={m:dict(type='endpoint',seeds=['20260916','20260917','20260918'],identity='frozen controlled family OOF') for m in ['baseline','copra']});c['methods']['zero']=dict(type='zero',seeds=['fixed'],identity='zero-change reference');(out/'config.json').write_text(json.dumps(c,indent=2));r=evaluate.run(out/'config.json');(out/'replay_summary.json').write_text(json.dumps(r,indent=2,allow_nan=False));assert (r['absolute']['common_units'],r['coverage']['common_pairs'],r['coverage']['local_pairs'])==(752,338,220)
print(json.dumps({'status':'PASS','units':len(u),'pairs':len(p),'source_links_used':sum(len(v) for v in links.values())}))
