import csv,json,tempfile,unittest,shutil
from pathlib import Path
from evaluate import run
SRC=Path(__file__).resolve().parent/'synthetic'
class Regression(unittest.TestCase):
 def setUp(self):
  self.t=tempfile.TemporaryDirectory(); self.p=Path(self.t.name)
  for n in ['units.csv','pairs.csv','predictions.csv','config.json']: shutil.copy(SRC/n,self.p/n)
 def tearDown(self): self.t.cleanup()
 def config(self): return json.loads((self.p/'config.json').read_text())
 def save(self,c): (self.p/'config.json').write_text(json.dumps(c))
 def test_disjoint(self):
  v=run(self.p/'config.json')['scopes']['local']; self.assertEqual(v['rank_common_landscapes'],0); self.assertNotIn('position_mean_spearman',v['contrasts']['B'])
 def test_partial(self):
  p=self.p/'predictions.csv'
  with p.open(newline='') as f: rows=list(csv.DictReader(f)); fields=list(rows[0])
  for r in rows:
   if r['method']=='B' and r['member_id'].startswith('L1_'): r['value']=-float(r['member_id'].split('P')[1])
  with p.open('w',newline='') as f: w=csv.DictWriter(f,fields); w.writeheader(); w.writerows(rows)
  v=run(self.p/'config.json')['scopes']['local']; c=v['contrasts']['B']['position_mean_spearman']; self.assertEqual((v['rank_common_landscapes'],v['rank_common_families']),(1,1)); self.assertAlmostEqual(c['reference_minus_method'],2); self.assertAlmostEqual(c['ci95'][0],2); self.assertEqual(c['families'],1)
 def test_multiplier(self):
  c=self.config(); c['methods']['A']['multiplier']='nan'; self.save(c)
  with self.assertRaisesRegex(ValueError,'multiplier'): run(self.p/'config.json')
  c['methods']['A']['multiplier']=1e308; self.save(c); r=run(self.p/'config.json'); self.assertLess(r['coverage']['successful_by_method']['A'],10)
 def test_endpoint_subtraction_overflow(self):
  c=self.config(); c['methods']={'A':dict(type='endpoint',seeds=['fixed'],identity='x')}; self.save(c)
  rows=[]
  with (self.p/'units.csv').open(newline='') as f:
   for u in csv.DictReader(f): rows.append(dict(method='A',seed='fixed',member_id=u['unit_id'],status='success',value=(-1e308 if u['unit_id'].endswith('WT') else 1e308),checkpoint='same'))
  with (self.p/'predictions.csv').open('w',newline='') as f: w=csv.DictWriter(f,list(rows[0])); w.writeheader(); w.writerows(rows)
  self.assertEqual(run(self.p/'config.json')['coverage']['successful_by_method']['A'],0)
if __name__=='__main__': unittest.main()
