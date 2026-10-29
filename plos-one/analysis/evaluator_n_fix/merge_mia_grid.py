"""Merge the N-safe MIA rerun into the main-grid MIA table.

Input:  ../mia_followup/results/mia_grid_with_followup.csv (v1.3.0 grid plus
        3,600-s follow-up; outcome of record), results/grid.csv (rerun of
        all 720 jobs with stored sequences, N-safe scoring, 9 in parallel), and
        results/grid_crash_retry.csv (every job that failed in grid.csv, rerun
        once, one job at a time; MIA crashes are intermittent, so a retry that
        returned a sequence replaces the failure, as in the 3,600-s follow-up).
Output: results/mia_grid_final.csv with the columns of the follow-up table;
        metrics come from the rerun. Outcome changes between the two runs
        are listed, and the rerun outcome is used, because only its
        sequences are retained.
"""

from __future__ import annotations

import csv
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
FOLLOWUP = os.path.join(HERE, "..", "mia_followup", "results", "mia_grid_with_followup.csv")
RERUN = os.path.join(HERE, "results", "grid.csv")
RETRY = os.path.join(HERE, "results", "grid_crash_retry.csv")
OUT = os.path.join(HERE, "results", "mia_grid_final.csv")
METRICS = ("theta_len", "n_ambiguous", "identity_pct", "recall", "precision", "f1")


def main() -> None:
    key = lambda r, d: (r["reference"], r["target"], int(r[d]), int(r["replicate"]))
    base = list(csv.DictReader(open(FOLLOWUP, newline="")))
    rerun = {key(r, "condition"): r for r in csv.DictReader(open(RERUN, newline="")) if r["method"] == "MIA"}
    retried = Counter()
    for r in csv.DictReader(open(RETRY, newline="")):
        retried[r["status"]] += 1
        if r["status"] == "completed":
            rerun[key(r, "condition")] = {**r, "retried": "1"}
    missing = [key(r, "depth") for r in base if key(r, "depth") not in rerun]
    if missing:
        raise SystemExit(f"rerun incomplete: {len(missing)} of {len(base)} jobs missing")
    changes = Counter()
    out = []
    for r in base:
        n = rerun[key(r, "depth")]
        status = "completed" if n["status"] == "completed" else "process_error"
        if status != r["final_status"]:
            changes[(r["final_status"], status)] += 1
        row = {**r, "final_status": status, "source": "nsafe_retry" if n.get("retried") else "nsafe_rerun",
               "final_runtime_s": n.get("runtime_s", ""),
               "theta_fasta": n.get("theta_fasta", ""), "theta_sha256": n.get("theta_sha256", "")}
        for m in METRICS:
            row[m] = n.get(m, "") if status == "completed" else ""
        out.append(row)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print("sequential retry of failed jobs:", dict(retried))
    print("final status:", Counter(r["final_status"] for r in out))
    print("outcome changes (follow-up -> rerun):", dict(changes) or "none")


if __name__ == "__main__":
    main()
