"""Compare a replay against the frozen full-scope reference values."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import fitz

ROOT = Path(__file__).resolve().parent
REF = json.loads((ROOT / "frozen/reference.json").read_text(encoding="utf-8"))
TOL = 1e-9


def read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


class Checker:
    def __init__(self):
        self.rows = []

    def add(self, name, ok, detail="", **extra):
        self.rows.append({"check": name, "status": "PASS" if ok else "FAIL",
                          "detail": detail, **extra})
        return ok

    def equal(self, name, actual, expected, tol=TOL):
        if isinstance(expected, dict):
            if not isinstance(actual, dict):
                return self.add(name, False, "type mismatch")
            ok = set(actual) == set(expected)
            worst = 0.0
            for key in expected:
                if key not in actual:
                    continue
                if isinstance(expected[key], (dict, list)):
                    sub = self.equal(f"{name}/{key}", actual[key], expected[key], tol)
                    ok = ok and sub
                elif isinstance(expected[key], (int, float)) and not isinstance(expected[key], bool):
                    if not isinstance(actual[key], (int, float)):
                        ok = False
                        continue
                    diff = abs(float(actual[key]) - float(expected[key]))
                    worst = max(worst, diff)
                    ok = ok and diff <= tol
                else:
                    ok = ok and actual[key] == expected[key]
            return self.add(name, ok, f"max numeric difference {worst:.3g}")
        if isinstance(expected, list):
            if not isinstance(actual, list) or len(actual) != len(expected):
                return self.add(name, False, "list length mismatch")
            ok = True
            worst = 0.0
            for index, (a, e) in enumerate(zip(actual, expected)):
                if isinstance(e, (dict, list)):
                    ok = self.equal(f"{name}[{index}]", a, e, tol) and ok
                elif isinstance(e, (int, float)) and not isinstance(e, bool):
                    diff = abs(float(a) - float(e))
                    worst = max(worst, diff)
                    ok = ok and diff <= tol
                else:
                    ok = ok and a == e
            return self.add(name, ok, f"max numeric difference {worst:.3g}")
        if isinstance(expected, (int, float)) and not isinstance(expected, bool):
            diff = abs(float(actual) - float(expected))
            return self.add(name, diff <= tol, f"difference {diff:.3g}")
        return self.add(name, actual == expected, "exact match")


def figure_record(path: Path) -> dict:
    with fitz.open(path) as doc:
        page = doc[0]
        words = sorted(w[4] for w in page.get_text("words"))
        return {
            "pages": doc.page_count,
            "page_size_pt": [round(page.rect.width, 4), round(page.rect.height, 4)],
            "text_tokens": len(words),
            "text_token_sha256": hashlib.sha256("\n".join(words).encode()).hexdigest(),
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", required=True, type=Path)
    parser.add_argument("--figures", type=Path)
    parser.add_argument("--external", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    replay = args.replay.resolve()
    check = Checker()

    primary = read(replay / "primary/replay_summary.json")
    check.add("primary/coverage", primary["coverage"]["common_pairs"] == 338
              and primary["absolute"]["common_units"] == 752
              and primary["coverage"]["local_pairs"] == 220,
              f"{primary['absolute']['common_units']} units / "
              f"{primary['coverage']['common_pairs']} pairs / "
              f"{primary['coverage']['local_pairs']} local")
    check.equal("primary/absolute", primary["absolute"]["methods"], REF["primary"]["absolute"]["methods"])

    table2 = read(replay / "checks/table2_keyed_check.json")
    check.add("table2/status", all(row["status"] == "PASS" for row in table2),
              f"{len(table2)} rows / {sum(len(r['actual']) for r in table2)} cells")
    check.equal("table2/values", [row["actual"] for row in table2],
                [row["actual"] for row in REF["table2_check"]])

    table3 = read(replay / "checks/table3_keyed_check.json")
    check.add("table3/status", all(row["status"] == "PASS" for row in table3),
              f"{len(table3)} rows")
    check.equal("table3/values", [row["actual"] for row in table3],
                [row["actual"] for row in REF["table3_check"]])

    methods = read(replay / "methods/method_results.json")
    ref_methods = REF["method_results"]
    for scope in ref_methods:
        check.equal(f"methods/{scope}/counts", {k: methods[scope].get(k) for k in ("units", "pairs", "bshape_pairs", "rank_pairs", "rank_landscapes")},
                    {k: ref_methods[scope].get(k) for k in ("units", "pairs", "bshape_pairs", "rank_pairs", "rank_landscapes")})
        check.equal(f"methods/{scope}/values", methods[scope]["methods"], ref_methods[scope]["methods"])

    pra = read(replay / "pra/verification.json")
    check.add("pra/all_pass", pra["all_pass"] is True, f"endpoints={pra['endpoints']}")
    check.add("pra/endpoints", pra["endpoints"] == REF["pra"]["endpoints"], str(pra["endpoints"]))
    check.add("pra/formula_error", pra["max_formula_error_G"] <= 1e-10, f"{pra['max_formula_error_G']:.3g}")
    check.add("pra/contribution_sum", pra["max_contribution_sum_error_G"] <= 1e-10,
              f"{pra['max_contribution_sum_error_G']:.3g}")

    m2 = read(replay / "mechanisms/m2_mechanism_check.json")
    check.add("m2/status", m2["status"] == "PASS" and m2["metrics"] == 72,
              f"metrics={m2['metrics']} max={m2['max_difference']:.3g}")
    three = read(replay / "mechanisms/three_arm_check.json")
    check.add("three_arm/status", three["status"] == "PASS" and three["cells"] == 30
              and three["seed_values"] == 90,
              f"cells={three['cells']} seed_values={three['seed_values']}")

    if args.figures:
        for number in range(1, 8):
            path = args.figures / f"figure{number}.pdf"
            if not path.is_file():
                check.add(f"figure{number}/exists", False, str(path))
                continue
            actual = figure_record(path)
            expected = REF["figures"][f"figure{number}"]
            ok = (actual["pages"] == expected["pages"]
                  and actual["page_size_pt"] == expected["page_size_pt"]
                  and actual["text_token_sha256"] == expected["text_token_sha256"])
            check.add(f"figure{number}/layout", ok,
                      f"actual={actual['page_size_pt']} tokens={actual['text_tokens']}")

    if args.external:
        external = read(args.external / "external_replay.json")
        check.add("external/status", external.get("status") == "PASS"
                  and external.get("cells_checked") == 64
                  and external.get("contrasts_checked") == 5,
                  f"cells={external.get('cells_checked')} contrasts={external.get('contrasts_checked')}")
        check.equal("external/family_metrics", external["family_metrics"],
                    REF["external"]["family_metrics"])
        check.equal("external/macro", external["macro"], REF["external"]["macro"])
        check.equal("external/contrasts", external["contrasts"], REF["external"]["contrasts"])

    failed = [row for row in check.rows if row["status"] == "FAIL"]
    report = {
        "status": "PASS" if not failed else "FAIL",
        "checks": check.rows,
        "failed": failed,
        "reference_version": REF["version"],
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "verification.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "checks": len(check.rows), "failed": len(failed)},
                     indent=2))
    if failed:
        for row in failed[:12]:
            print(json.dumps(row))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
