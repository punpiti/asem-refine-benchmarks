"""Compare the tie-rule reruns with the released asem-refine-benchmarks v1.3.0 results.

The reruns (rerun_tie_rule.sh / rerun_real_wgs.sh in ~/plos-reruns/bench-v1.3.0)
regenerated result sets that predate, or cannot be shown to postdate, the
2026-09-30 tie rule (tied insertions omitted, ties keep the current base),
using an unmodified git-archive export of v1.3.0. Rows are matched on every
non-metric key column; reported per set: matched rows, rows whose F1 changed,
largest |delta F1| and |delta identity|. The rerun CSVs are copied to results/.
"""

from __future__ import annotations

import csv
import os
import shutil

BENCH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
HERE = os.path.dirname(os.path.abspath(__file__))
SETS = {
    "read_length_experiment": "results.csv",
    "parameter_sensitivity": "results.csv",
    "circular_boundary": "results.csv",
    "real_wgs_benchmark": "grid_results.csv",
    "real_wgs_hybrid": "final_results.csv",
    "real_wgs_hybrid_boundary": "final_results.csv",
}
METRICS = {"identity_pct", "recall", "precision", "f1", "c_theta", "w_theta", "u_t", "runtime_s", "wall_s",
           "theta_len", "n_contigs_total", "n_anchors_total", "last_iteration", "hit_iteration_cap",
           "host", "error", "status"}


def main() -> None:
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    lines = []
    for name, fname in SETS.items():
        old_rows = list(csv.DictReader(open(os.path.join(BENCH, "released", name, fname), newline="")))
        new_path = os.path.join(BENCH, "scripts", "results", name, fname)
        new_rows = list(csv.DictReader(open(new_path, newline="")))
        shutil.copy(new_path, os.path.join(HERE, "results", f"{name}__{fname}"))
        keys = [c for c in old_rows[0] if c not in METRICS and not any(m in c for m in ("f1", "identity", "recall", "precision", "theta", "_s"))]
        key = lambda r: tuple(r[c] for c in keys)
        old = {key(r): r for r in old_rows}
        new = {key(r): r for r in new_rows}
        common_keys = old.keys() & new.keys()
        f1_cols = [c for c in old_rows[0] if c.endswith("f1")]
        id_cols = [c for c in old_rows[0] if c.endswith("identity_pct")]

        def delta(cols):
            out = []
            for k in common_keys:
                for c in cols:
                    if old[k][c] not in ("", None) and new[k].get(c) not in ("", None):
                        out.append(abs(float(new[k][c]) - float(old[k][c])))
            return out
        df, di = delta(f1_cols), delta(id_cols)
        changed = sum(1 for k in common_keys for c in f1_cols
                      if old[k][c] and new[k].get(c) and float(old[k][c]) != float(new[k][c]))
        lines.append(f"{name}: released {len(old)} rows, rerun {len(new)}, matched {len(common_keys)} "
                     f"(keys: {','.join(keys)}); F1 changed in {changed}; "
                     f"max |dF1| {max(df, default=0):.2g}; max |dIdentity| {max(di, default=0):.2g} points")
    text = "\n".join(lines) + "\n"
    open(os.path.join(HERE, "results", "comparison.txt"), "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
