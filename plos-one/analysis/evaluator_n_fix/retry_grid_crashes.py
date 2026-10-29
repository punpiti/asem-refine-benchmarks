"""Retry, one job at a time, every main-grid MIA job that crashed in the N-safe rerun.

The N-safe grid rerun ran 9 jobs in parallel. Two jobs that completed in the
v1.3.0 grid crashed there (MIA exit -6/-11), so this checks whether the MIA
crashes are reproducible on an otherwise idle machine.

Writes results/grid_crash_retry.csv (same columns as results/grid.csv).
"""

from __future__ import annotations

import csv
import os

from rerun_external import FIELDS, RESULTS_DIR, grid_job


def main() -> None:
    with open(os.path.join(RESULTS_DIR, "grid.csv"), newline="") as fh:
        failed = [r for r in csv.DictReader(fh) if r["status"] != "completed"]
    out_csv = os.path.join(RESULTS_DIR, "grid_crash_retry.csv")
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        for n, r in enumerate(failed, start=1):
            job = ("grid", r["method"], r["reference"], r["target"], int(r["condition"]),
                   int(r["replicate"]), r["n_reads"])
            row = grid_job(job)
            writer.writerow({k: row.get(k, "") for k in FIELDS})
            fh.flush()
            print(f"[{n}/{len(failed)}] {job[2]} {job[3]} {job[4]}X r{job[5]} {row['status']} "
                  f"runtime={row.get('runtime_s')} f1={row.get('f1', '')}", flush=True)


if __name__ == "__main__":
    main()
