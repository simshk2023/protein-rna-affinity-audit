import sys,unittest
from pathlib import Path
p=Path(__file__).resolve().parent/'toolkit'
sys.path.insert(0,str(p))
suite=unittest.defaultTestLoader.discover(str(p/'tests'))
raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
