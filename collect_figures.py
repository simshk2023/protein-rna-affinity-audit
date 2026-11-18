"""Copy the archived v5 renderers' own page output as the manuscript figures.

The current v5 manuscript uses the renderer pages directly.  The optional
composition helper is kept for the supplementary S14 layout only.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

MAPPING = {
    1: ("figure1", "rendered/variants/glossy/molecular_context_grounded.pdf"),
    2: ("figure2", "figure1.pdf"),
    3: ("figure3", "figure2.pdf"),
    4: ("figure4", "figure3.pdf"),
    5: ("figure5", "figure4.pdf"),
    6: ("figure6", "figure5.pdf"),
    7: ("figure7", "figure6.pdf"),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    copied = {}
    for number, (directory, name) in MAPPING.items():
        source = args.inputs / directory / name
        if not source.is_file():
            raise ValueError(f"missing rendered Figure {number}: {source}")
        target = args.output / f"figure{number}.pdf"
        shutil.copyfile(source, target)
        copied[number] = {"path": str(target), "bytes": target.stat().st_size}
    (args.output / "figure_manifest.json").write_text(json.dumps(copied, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(copied, indent=1))


if __name__ == "__main__":
    main()
