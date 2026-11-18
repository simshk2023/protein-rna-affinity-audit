"""Recommended offline route. All real outputs stay in the user-selected restricted directory."""
import argparse,subprocess,sys,json,time
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--inputs',required=True);ap.add_argument('--output',required=True);ap.add_argument('--skip-figures',action='store_true');a=ap.parse_args();h=Path(__file__).resolve().parent;I=Path(a.inputs).resolve();O=Path(a.output).resolve();O.mkdir(parents=True,exist_ok=True);log=[]
def run(script,args):
 t=time.time();q=subprocess.run([sys.executable,str(h/script),*map(str,args)],capture_output=True,text=True,encoding='utf-8',errors='replace');(O/(script+'.log')).write_text(q.stdout+q.stderr,encoding='utf-8');log.append(dict(entry=script,exit_code=q.returncode,seconds=time.time()-t));(O/'run_manifest.json').write_text(json.dumps(log,indent=2));print(log[-1],flush=True)
 if q.returncode:raise SystemExit(q.stderr[-1500:])
run('import_authorized.py',['--pronab',I/'pronab.csv','--predictions',I/'predictions.csv','--output',O/'primary'])
run('check_table2.py',[O/'primary/replay_summary.json','--manuscript',h/'main_tables_snapshot.tex','--output',O/'checks'])
if (I/'services/service_pair_predictions.json').is_file():
 shutil=__import__('shutil');(O/'services').mkdir(exist_ok=True);shutil.copyfile(I/'services/service_pair_predictions.json',O/'services/service_pair_predictions.json')
else:
 run('replay_services.py',['--inputs',I/'services','--output',O/'services'])
run('replay_methods.py',['--inputs',I/'methods','--service-predictions',O/'services/service_pair_predictions.json','--output',O/'methods'])
run('check_table3.py',[O/'methods/method_results.json','--manuscript',h/'main_tables_snapshot.tex','--output',O/'checks/table3_keyed_check.json'])
run('replay_pra_mechanism.py',['--inputs',I/'mechanisms','--output',O/'pra','--formula',I/'pra_formula.json'])
run('replay_mechanisms.py',['--inputs',I/'mechanisms','--methods',I/'methods','--figures',I/'figures','--output',O/'mechanisms'])
if not a.skip_figures:run('build_figures.py',['--inputs',I/'figures','--output',O/'figures'])
print('Offline delivery route PASS',flush=True)
