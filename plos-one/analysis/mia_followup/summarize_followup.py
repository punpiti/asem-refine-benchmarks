"""Acceptance check and merged view of the MIA follow-up.

Reads the released v1.3.0 MIA grid and followup_results.csv and writes:

  results/acceptance.json            expected/observed counts, duplicate and
                                     missing-key checks, status counts,
                                     runtime distributions, tool versions
  results/mia_grid_with_followup.csv one row per original job key: the
                                     original row's status plus the outcome of
                                     its largest-budget follow-up attempt

Completed accuracy, censoring and process errors stay in separate columns; a
job without scoreable output is never scored as F1 = 0.
"""

from __future__ import annotations

import csv
import glob
import json
import os
import statistics
from collections import Counter

import run_mia_followup as runner

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")


def _quantiles(values: list[float]) -> dict:
    if not values:
        return {"n": 0}
    values = sorted(values)
    q = statistics.quantiles(values, n=4) if len(values) > 1 else [values[0]] * 3
    return {"n": len(values), "min": round(values[0], 2), "q1": round(q[0], 2),
            "median": round(q[1], 2), "q3": round(q[2], 2), "max": round(values[-1], 2)}


def main() -> None:
    original = runner.load_original()
    original_keys = [runner._key(r) for r in original]
    with open(runner.FOLLOWUP_CSV, newline="") as fh:
        followup = list(csv.DictReader(fh))

    attempt_keys = [(*runner._key(r), int(r["timeout_budget_s"])) for r in followup]
    latest: dict[tuple, dict] = {}
    for row in followup:
        key = runner._key(row)
        if key not in latest or int(row["timeout_budget_s"]) > int(latest[key]["timeout_budget_s"]):
            latest[key] = row

    expected = {runner._key(r) for r in original if r["original_status"] in {"timeout", "process_error"}}
    missing = sorted(expected - set(latest))
    unexpected = sorted(set(latest) - expected)

    merged = []
    for row in original:
        key = runner._key(row)
        follow = latest.get(key)
        final = follow if follow else row
        final_status = follow["status"] if follow else row["original_status"]
        merged.append({
            "reference": row["reference"], "target": row["target"],
            "divergence": row["divergence"], "depth": row["depth"], "replicate": row["replicate"],
            "original_status": row["original_status"],
            "original_runtime_s": row["runtime_s"],
            "followup_budget_s": follow["timeout_budget_s"] if follow else "",
            "final_status": final_status,
            "final_runtime_s": follow["wall_s"] if follow else row["runtime_s"],
            "source": "followup" if follow else "v1.3.0",
            **{m: (final[m] if final_status == "completed" else "")
               for m in ("theta_len", "n_ambiguous", "identity_pct", "recall", "precision", "f1")},
        })
    with open(os.path.join(RESULTS, "mia_grid_with_followup.csv"), "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(merged[0]))
        writer.writeheader()
        writer.writerows(merged)

    manifests = []
    for path in sorted(glob.glob(os.path.join(runner.MANIFEST_DIR, "*.json"))):
        with open(path) as fh:
            m = json.load(fh)
        manifests.append({k: m[k] for k in ("selection", "timeout_s", "parallel", "n_jobs", "host",
                                             "created_utc", "mia_sha256", "ma_sha256", "mia_package",
                                             "benchmark_repo_commit", "benchmark_repo_dirty_scripts",
                                             "original_csv_sha256")}
                         | {"file": os.path.relpath(path, HERE)})

    followup_completed = [float(r["wall_s"]) for r in latest.values() if r["status"] == "completed"]
    acceptance = {
        "original_rows": len(original),
        "original_unique_keys": len(set(original_keys)),
        "original_status_counts": dict(Counter(r["original_status"] for r in original)),
        "expected_followup_keys": len(expected),
        "followup_attempts": len(followup),
        "followup_duplicate_attempts": len(attempt_keys) - len(set(attempt_keys)),
        "followup_keys_covered": len(set(latest) & expected),
        "missing_followup_keys": [list(k) for k in missing],
        "unexpected_followup_keys": [list(k) for k in unexpected],
        "input_mismatches": sum(r["status"] == "input_mismatch" for r in followup),
        "followup_status_by_original_status": {
            status: dict(Counter(latest[runner._key(r)]["status"] for r in original
                                 if r["original_status"] == status and runner._key(r) in latest))
            for status in ("timeout", "process_error")
        },
        "final_status_counts": dict(Counter(r["final_status"] for r in merged)),
        "runtime_s": {
            "original_completed": _quantiles([float(r["runtime_s"]) for r in original
                                              if r["original_status"] == "completed"]),
            "followup_completed": _quantiles(followup_completed),
            "followup_peak_rss_kb": _quantiles([float(r["peak_rss_kb"]) for r in latest.values()
                                                if r["peak_rss_kb"]]),
        },
        "manifests": manifests,
        "accepted": not missing and not unexpected
                    and len(set(original_keys)) == len(original) == 720
                    and len(attempt_keys) == len(set(attempt_keys))
                    and not any(r["status"] == "input_mismatch" for r in followup),
    }
    with open(os.path.join(RESULTS, "acceptance.json"), "w") as fh:
        json.dump(acceptance, fh, indent=2)
    print(json.dumps({k: v for k, v in acceptance.items() if k != "manifests"}, indent=2))


if __name__ == "__main__":
    main()
