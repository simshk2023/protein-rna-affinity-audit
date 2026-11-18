"""Rebuild the v5 display layout from the archived figure renderers."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import fitz


# Rectangles are in PDF points, measured down and right from the top-left.
COMPOSITIONS = {
    "figure3.pdf": (
        (506, 268),
        [
            ("3", (0, 0, 192, 265), (0, 1, 192, 266)),
            ("3", (0, 417.5, 198, 432), (207, 1, 405, 15.5)),
            ("3", (14, 438, 251, 560), (216, 20, 453, 142)),
            ("3", (312, 438, 498, 560), (258, 143, 444, 265)),
        ],
    ),
    "figure6.pdf": (
        (244, 416.17728),
        [
            ("6", (8, 0, 27, 199), (7, 0, 24.3413, 181.6267)),
            ("6", (27, 0, 260, 210), (24.3413, 0, 237, 191.667)),
            ("alignment", (0, 104, 243.77953, 320.31497),
             (0, 199.667, 244, 416.17728)),
        ],
    ),
    "figureS14.pdf": (
        (504.57, 211),
        [
            ("6", (0, 0, 267, 201.5), (0, 0, 267, 201.5)),
            ("6", (35, 201.5, 259, 210), (35, 201.5, 259, 210)),
            ("6", (267, 0, 504.566925, 201.5), (267.00156, 0, 504.56848, 201.5)),
        ],
    ),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archived_renders", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--figures", default="1,2,3,4,5,6,7,S14")
    args = parser.parse_args()
    figures = {item.strip().upper() for item in args.figures.split(",")}
    if not figures <= {"1", "2", "3", "4", "5", "6", "7", "S14"}:
        parser.error("--figures accepts a comma-separated list from 1–7 and S14")

    args.output.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parent
    assets = {}
    for figure in json.loads((root / "figure_manifest.json").read_text(encoding="utf-8")):
        number = str(figure["figure"])
        assets[number] = args.archived_renders / f"figure{number}" / figure["expected_pdf"]
        if number in figures and number not in ("3", "6"):
            shutil.copyfile(assets[number], args.output / f"figure{number}.pdf")

    if "6" in figures:
        work = args.output / "alignment_source"
        work.mkdir(exist_ok=True)
        for name in ("aggregates.json", "draw_alignment.py"):
            shutil.copyfile(root / "manuscript_v5" / name, work / name)
        subprocess.run([sys.executable, "-B", str(work / "draw_alignment.py")], check=True)
        assets["alignment"] = work / "figure6b_source.pdf"

    for name, ((width, height), panels) in COMPOSITIONS.items():
        if name.removeprefix("figure").removesuffix(".pdf") not in figures:
            continue
        with fitz.open() as output:
            page = output.new_page(width=width, height=height)
            for key, clip, destination in panels:
                with fitz.open(assets[key]) as source:
                    page.show_pdf_page(
                        fitz.Rect(destination), source, 0, clip=fitz.Rect(clip),
                        keep_proportion=False,
                    )
            output.save(args.output / name, garbage=4, deflate=True)
    print(f"Composed figures {args.figures}.")


if __name__ == "__main__":
    main()
