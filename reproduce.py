"""Reproduce the CoPRA v5 results from the public package plus reader imports.

Modes
  smoke      run the bundled synthetic tests only
  prepare    build work/inputs from the package and the imported third-party files
  primary    Table 2 + Supplementary Table S36 (752 units / 338 pairs / 220 local)
  full       Table 2, Table 3, S36, S37, PRA, M2, three arms and Figures 1-7/S14
  external   Supplementary S13 and Figure E1 from imported observations
  all        full + external + the frozen-value verification
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run(script: str, *args, cwd: Path | None = None) -> None:
    subprocess.run([sys.executable, "-B", str(ROOT / script), *map(str, args)],
                   check=True, cwd=cwd or ROOT)


def verify_package() -> int:
    manifest = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
    checked = 0
    for relative, record in manifest["files"].items():
        path = ROOT / relative
        if not path.is_file():
            raise ValueError(f"missing package file: {relative}")
        with path.open("rb") as fh:
            digest = hashlib.file_digest(fh, "sha256").hexdigest()
        if digest != record["sha256"] or path.stat().st_size != record["bytes"]:
            raise ValueError(f"package file differs from manifest: {relative}")
        checked += 1
    return checked


def prepare(args) -> Path:
    work = (args.work or (args.output / "work")).resolve()
    cmd = [sys.executable, "-B", str(ROOT / "prepare_inputs.py"),
           "--pronab", str(args.pronab.resolve()),
           "--output", str(work)]
    if args.pra_formula:
        cmd += ["--pra-formula", str(args.pra_formula.resolve())]
    if args.structures:
        cmd += ["--structures", str(args.structures.resolve())]
    subprocess.run(cmd, check=True, cwd=ROOT)
    return work


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("smoke", "prepare", "primary", "full", "external", "all"))
    parser.add_argument("--pronab", type=Path)
    parser.add_argument("--pra-formula", type=Path)
    parser.add_argument("--structures", type=Path)
    parser.add_argument("--observations", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--work", type=Path)
    args = parser.parse_args()

    output = args.output.resolve()
    if output.is_relative_to(ROOT):
        parser.error("--output must be outside the release source directory")
    if args.mode == "smoke":
        run("run_tests.py")
        run("toolkit/evaluate.py", "toolkit/example/config.json", "--output", output)
        return

    if args.mode == "external" and not args.observations:
        parser.error("--observations is required for the external route")
    if args.mode in ("primary", "full", "all") and not args.pronab:
        parser.error("--pronab is required; see study_outputs/primary/ACCESS.md")
    if args.mode in ("full", "all") and not args.pra_formula:
        parser.error("--pra-formula is required for Figure 5 and PRA verification")

    checked = verify_package()
    print(f"Package manifest verified: {checked} files.", flush=True)

    if args.mode == "external":
        output.mkdir(parents=True, exist_ok=True)
        run("replay_open_external.py", "--observations", args.observations.resolve(),
            "--output", output)
        return

    work = prepare(args) if args.work is None else args.work.resolve()
    output.mkdir(parents=True, exist_ok=True)

    if args.mode in ("primary", "prepare"):
        if args.mode == "prepare":
            return
        run("import_authorized.py", "--pronab", work / "pronab.csv",
            "--predictions", ROOT / "study_outputs/primary/predictions.csv",
            "--output", output / "primary")
        run("check_table2.py", output / "primary/replay_summary.json",
            "--manuscript", ROOT / "main_tables_snapshot.tex", "--output", output / "checks")
        sys.path.insert(0, str(ROOT))
        from check_v5_tables import check_primary
        check_primary(output)
        print("Main Table 2 and Supplementary Table S36 PASS.", flush=True)
        return

    full = output / "full"
    run("run_delivery.py", "--inputs", work, "--output", full)
    run("collect_figures.py", "--inputs", full / "figures",
        "--output", output / "manuscript_figures")
    run("compose_v5.py", full / "figures", "--output", output / "manuscript_figures",
        "--figures", "S14")
    run("verify_frozen.py", "--replay", full, "--figures", output / "manuscript_figures",
        "--output", output / "verification")
    if args.observations:
        run("replay_open_external.py", "--observations", args.observations.resolve(),
            "--output", output / "external")
        run("verify_frozen.py", "--replay", full, "--figures", output / "manuscript_figures",
            "--external", output / "external", "--output", output / "verification")
    print(f"Full replay and frozen-value verification PASS: {output}", flush=True)


if __name__ == "__main__":
    main()
