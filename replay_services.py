"""Offline archived service responses to endpoint/pair mappings. No network imports."""
import argparse,json,csv,math
from pathlib import Path
from collections import defaultdict
from decimal import Decimal
ap=argparse.ArgumentParser();ap.add_argument('--inputs',required=True);ap.add_argument('--output',required=True);a=ap.parse_args();I=Path(a.inputs);O=Path(a.output);O.mkdir(parents=True,exist_ok=True)
def rd(p,delimiter=','):return list(csv.DictReader(p.open(encoding='utf-8-sig',newline=''),delimiter=delimiter))
def body(p):
 b=p.read_bytes()
 try:return json.loads(b.decode('utf-8-sig'))
 except UnicodeDecodeError:return json.loads(b.decode('cp1252'))
P=I/'panther';scores={}
for x in json.loads((P/'ledger.json').read_text()):
 j=body(P/x['body']);assert j['success'] is True
 assert j['uploaded_file']['sha256'].lower()==x['input_sha256'].lower()
 scores[x['request_id']]=Decimal(str(j['results'][0]['score']))
pan={r['pair_id']:float(-(scores[r['mut_request_id']]-scores[r['wt_request_id']])) for r in rd(P/'cohort_pairs.csv')};assert len(scores)==110 and len(pan)==87
M=I/'mcsm';result=defaultdict(list)
for req in rd(M/'request_manifest.csv'):
 for r in rd(M/'responses'/(req['request_id']+'.tsv'),delimiter='\t'):
  result[(r['PDB_FILE'].upper(),r['CHAIN'],r['WILD_RES'],int(r['RES_POS']),r['MUT_RES'])].append(float(r['PRED_DDG']))
mc={}
for r in rd(M/'scope_348.csv'):
 if r['pdb_id_applicable']!='True':continue
 k=(r['pdb'].upper(),r['chain'],r['wt1'],int(r['pdb_resseq']),r['mut1']);v=result.get(k,[])
 if len(v)!=1:raise ValueError('response join is not unique: '+r['pair_id'])
 mc[r['pair_id']]=v[0]/1.3635606572452885
assert len(mc)==267
(O/'service_pair_predictions.json').write_text(json.dumps(dict(panther=pan,mcsm=mc),indent=2));print(json.dumps({'status':'PASS','panther_endpoint_responses':len(scores),'panther_pairs':len(pan),'mcsm_pairs':len(mc)}))
