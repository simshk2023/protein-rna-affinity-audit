import argparse,json,re
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('results');ap.add_argument('--manuscript',required=True);ap.add_argument('--output',required=True);a=ap.parse_args();r=json.loads(Path(a.results).read_text());tex=Path(a.manuscript).read_text(encoding='utf-8-sig');start=tex.index('Method or matched reference');panther=tex.index('PANTHER native-score', start);end=tex.find(r'\end{tabularx}', panther);b=tex[start:end if end != -1 else len(tex)];rows=[];scope=None
aliases={'CoPRA (controlled)':'copra','Additive':'baseline','CoPRA additive (OOF)':'baseline','DeePNAP':'deeppnap','Sequence additive':'seq_baseline','OmniBioTE official':'omni_official','OmniBioTE adapted':'omni_controlled','Zero change':'zero','PRdeltaGPred':'prdelta','mCSM-NA':'mcsm','GBNL reconstruction (OOF)':'gbnl','SVSBI reconstruction (OOF)':'svsbi','ProRB author checkpoint':'prorb','PRA-Pred':'pra','M2-PRNet, indices aligned':'m2_aligned','PANTHER':'panther'}
for line in b.splitlines():
 if 'Common scope' in line:scope=re.search(r'; (S\d+)',line)[1];continue
 if 'PANTHER native-score' in line:scope='S11';continue
 if '&' not in line or scope is None:continue
 cells=[x.strip().replace('\\\\','').strip() for x in line.split('&')];name=cells[0]
 if name not in aliases:continue
 method=aliases[name];x=r[scope]['methods'][method]
 if scope=='S11':
  expected=[float(re.search(r'\d+',cells[1])[0]),float(re.search(r'[\d.]+',cells[2])[0]),float(cells[3])];actual=[x['zeros'],x['strict_direction'],x['rank']];rounded=[round(actual[0]),round(actual[1],1),round(actual[2],4)]
 else:
  expected=[None if c=='---' else float(c) for c in cells[1:4]];actual=[x.get('absolute'),x['ordinary'],x['cmae']];rounded=[None if v is None else round(v,4) for v in actual]
 good=all(e is None or abs(e-v)<1e-10 for e,v in zip(expected,rounded));rows.append(dict(scope=scope,method=method,expected=expected,actual=actual,status='PASS' if good else 'FAIL'))
Path(a.output).write_text(json.dumps(rows,indent=2));assert len(rows)==37,len(rows)
assert all(x['status']=='PASS' for x in rows),[x for x in rows if x['status']!='PASS'];print(json.dumps({'status':'PASS','rows':len(rows)}))
