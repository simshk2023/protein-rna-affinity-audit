# Reusable affinity and mutation-effect evaluation

A lightweight, source-qualified evaluation of absolute affinity, paired mutation effects and within-landscape response. This bundle uses synthetic data only; the project-specific replay uses separately held restricted sources. It does not train or run a model.

## Quick start

Requires Python 3.10+ and NumPy. From any working directory, substitute the absolute path to this bundle:

```powershell
python C:/path/to/bundle/evaluate.py C:/path/to/bundle/example/config.json --output C:/path/to/output
```

The supplied example is entirely synthetic. For your own data, prepare the three CSVs and JSON configuration described in [SCHEMA.md](SCHEMA.md), run the same command with your configuration, and protect outputs containing family identities or restricted measurements. `summary.json` contains aggregate family-level scores and coverage; `report.md` is a compact overview. The core imports as `evaluate` and exposes `run(config_path)` and `landscape_metrics`.

Read [PROTOCOL.md](PROTOCOL.md) before interpreting or comparing results. [REPLAY.md](REPLAY.md) links the study-specific replay, and the package-level [README](../README.md) gives the smoke command.
