"""Stage authorized figure data and execute archived renderers in a separate output tree."""
import argparse,json,shutil,subprocess,sys,time
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--inputs',required=True);ap.add_argument('--output',required=True);ap.add_argument('--figures',default='1,2,3,4,5,6,7');a=ap.parse_args();h=Path(__file__).resolve().parent;o=Path(a.output).resolve();o.mkdir(parents=True,exist_ok=True);results=[]
sources_path=h/'study_outputs/figures/sources.json';sources=json.loads(sources_path.read_text()) if sources_path.is_file() else {}
for m in json.loads((h/'figure_manifest.json').read_text()):
 if str(m['figure']) not in a.figures.split(','):continue
 d=o/f"figure{m['figure']}";shutil.copytree(Path(a.inputs)/d.name,d,dirs_exist_ok=True);shutil.copytree(h/'figure_code'/d.name,d,dirs_exist_ok=True)
 t=time.time();args=[m['entry'],*sources.get(str(m['figure']),{}).get('renderer_args',[])];r=subprocess.run([sys.executable,*args],cwd=d,capture_output=True,text=True,encoding='utf-8',errors='replace');(d/'build.log').write_text(r.stdout+r.stderr,encoding='utf-8');results.append(dict(figure=m['figure'],exit_code=r.returncode,seconds=time.time()-t,pdf_exists=(d/m['expected_pdf']).exists()));(o/'build_results.json').write_text(json.dumps(results,indent=2));print(results[-1],flush=True)
 if r.returncode:print(r.stderr[-1500:],flush=True)
assert all(r['exit_code']==0 and r['pdf_exists'] for r in results),'figure build failure; see individual logs'
