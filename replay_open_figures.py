"""Render manuscript figures from shipped author layers and authorized local inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / "study_outputs/figures"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figures", default="1,3,4,6,S14")
    parser.add_argument("--inputs", type=Path, help="Authorized figures directory for Figures 2, 5 and 7")
    args = parser.parse_args()
    output = args.output.resolve()
    if output.is_relative_to(ROOT):
        parser.error("--output must be outside the package")
    sources = json.loads((SOURCES / "sources.json").read_text())
    requested = list(dict.fromkeys(x.strip().upper() for x in args.figures.split(",")))
    if not requested or any(x not in {*sources, "S14"} for x in requested):
        parser.error("--figures uses comma-separated values from 1–7 and S14")
    sources_needed = {"6" if value == "S14" else value for value in requested}
    frozen = {r["path"]: r for r in json.loads((ROOT / "input_manifest.json").read_text())}
    for number in sources_needed:
        for rel in sources[number]["authorized"]:
            if args.inputs is None or not (args.inputs / rel).is_file():
                parser.error(f"Figure {number} needs matching authorized figure input {rel}; see FIGURES.md")
            source = args.inputs / rel
            record = frozen["figures/" + rel]
            if source.stat().st_size != record["bytes"] or hashlib.sha256(source.read_bytes()).hexdigest() != record["sha256"]:
                parser.error(f"Figure {number} authorized input differs from the recorded snapshot: {rel}")
    entries = json.loads((ROOT / "figure_manifest.json").read_text())
    with tempfile.TemporaryDirectory(prefix="copra_figures_") as work:
        stages = Path(work)
        for entry in entries:
            number = str(entry["figure"])
            if number not in sources_needed:
                continue
            stage = stages / f"figure{number}"
            shutil.copytree(ROOT / "figure_code" / stage.name, stage)
            spec = sources[number]
            for key, origin in (("bundled", SOURCES), ("authorized", args.inputs)):
                for rel in spec[key]:
                    source = origin / rel
                    target = stages / rel
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(source, target)
            subprocess.run([sys.executable, "-B", entry["entry"], *spec.get("renderer_args", [])], cwd=stage, check=True)
        subprocess.run([sys.executable, "-B", str(ROOT / "compose_v5.py"), str(stages),
                        "--output", str(output), "--figures", ",".join(requested)], check=True)


if __name__ == "__main__":
    main()
