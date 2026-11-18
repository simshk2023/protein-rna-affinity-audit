import csv,json,tempfile,unittest
from pathlib import Path
from evaluate import run
class Integration(unittest.TestCase):
 def test_missing_and_score_only(self):
  source=Path(__file__).resolve().parents[1]/'example'
  with tempfile.TemporaryDirectory() as td:
   dest=Path(td)
   for name in ['units.csv','pairs.csv','predictions.csv','config.json']: (dest/name).write_bytes((source/name).read_bytes())
   c=json.loads((dest/'config.json').read_text()); c['methods']={'score':c['methods']['score']}; c['reference']='score'; (dest/'config.json').write_text(json.dumps(c))
   r=run(dest/'config.json'); self.assertEqual(r['coverage']['common_pairs'],10); metrics=r['scopes']['local']['methods']['score']['metrics']; self.assertNotIn('ordinary_mae',metrics); self.assertNotIn('cmae',metrics); self.assertEqual(r['absolute']['common_units'],0)
   with (dest/'predictions.csv').open(newline='',encoding='utf8') as f: rows=list(csv.DictReader(f)); fields=list(rows[0])
   row=next(x for x in rows if x['method']=='score'); row['value']='nan'
   with (dest/'predictions.csv').open('w',newline='',encoding='utf8') as f: w=csv.DictWriter(f,fields); w.writeheader(); w.writerows(rows)
   r=run(dest/'config.json'); self.assertEqual(r['coverage']['common_pairs'],9); self.assertEqual(r['coverage']['local_pairs'],5)
if __name__=='__main__': unittest.main()
