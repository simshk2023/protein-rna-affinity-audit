"""Fetch and verify the four wwPDB coordinate files used by Figure 1.

    python -B third_party/fetch_structures.py --output ../structures

Files are CC0 1.0.  The script refuses a file whose SHA-256 differs from the
recorded study input; an entry ID alone is not treated as sufficient proof of
byte identity.
"""
import argparse, hashlib, json, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CHECKS = json.loads((ROOT / "checksums.json").read_text())["structures"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--base", default="https://files.rcsb.org/download")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for pdb, record in CHECKS.items():
        target = args.output / f"{pdb}.pdb"
        url = f"{args.base}/{pdb}.pdb"
        print(f"fetch {url}")
        with urllib.request.urlopen(url) as response:
            data = response.read()
        digest = hashlib.sha256(data).hexdigest()
        if digest != record["sha256"] or len(data) != record["bytes"]:
            raise SystemExit(f"{pdb}: downloaded bytes differ from the recorded study input")
        target.write_bytes(data)
        print(f"  {target} sha256={digest}")


if __name__ == "__main__":
    main()
