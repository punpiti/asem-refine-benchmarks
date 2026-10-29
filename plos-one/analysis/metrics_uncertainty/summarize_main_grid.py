"""Paired, pair-clustered summaries of the main Phase-1 grid (Phase 4).

Grain and estimands (fixed before the summaries were inspected):
  * Job = (reference, target, depth, replicate). Methods share the reference,
    target, depth, replicate and simulated reads, so method comparisons are
    paired job differences (method minus ASEM).
  * Biological unit = species pair. Replicates are technical. Ordered pairs
    A->B and B->A share both genomes, so uncertainty resamples UNORDERED
    species pairs (clusters).
  * same-genus and same-family tiers have only 2 ordered pairs (1 unordered
    pair) each: report pair-level means and replicate ranges, no interval.
  * same-order tier (26 ordered pairs, 13 unordered clusters): descriptive
    cluster-bootstrap 95% percentile interval of the mean paired difference,
    10,000 resamples, seed 20261028. Each cluster contributes the mean of its
    jobs. No significance tests.
  * Accuracy is summarised over completed jobs only; output status is a
    separate table. A missing output is never scored as F1 = 0.
  * Internal methods: final iteration row. MIA: all 720 jobs rerun with the
    3,600-s limit, sequences stored and scored with the N-safe evaluator
    (plos-one/analysis/evaluator_n_fix/merge_mia_grid.py).
  * NOVOPlasty and GetOrganelle: scores after strand and circular-origin
    normalization (normalized_successful_rerun.csv), as in the main text and
    figures. Within 1-8X no stored NOVOPlasty/GetOrganelle output changes
    under the N-safe evaluator (evaluator_n_fix/rescore_stored_novoplasty.py).

Outputs in results/: status_by_cell.csv, accuracy_by_cell.csv,
paired_vs_asem.csv, pair_level.csv.
"""

from __future__ import annotations

import csv
import os
import random
import statistics as st
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
GRIDS = os.path.join(ROOT, "scripts", "results")
MIA_MERGED = os.path.join(ROOT, "plos-one", "analysis", "evaluator_n_fix", "results", "mia_grid_final.csv")
OUT = os.path.join(HERE, "results")

METHODS = {
    "ASEM": "baseline1_ieee_access",
    "ASEM_recursive": "baseline2_ojemb",
    "ASEM_Hybrid": "baseline_hybrid_boundary",
    "SRSC": "baseline3_ecticon",
    "Pilon": "ext_pilon",
    "MIA": None,
    "NOVOPlasty": "ext_novoplasty",
    "GetOrganelle": "ext_getorganelle",
}
DEPTHS = [str(d) for d in range(1, 9)]
N_BOOT = 10_000
SEED = 20261028


def _key(r: dict) -> tuple:
    return (r["reference"], r["target"], r["depth"], r["replicate"])


def load(method: str) -> dict[tuple, dict]:
    """Return {job key: {divergence, status, f1, identity}} for the 1-8X grid."""
    out: dict[tuple, dict] = {}
    if method == "MIA":
        for r in csv.DictReader(open(MIA_MERGED, newline="")):
            done = r["final_status"] == "completed"
            out[_key(r)] = {"divergence": r["divergence"], "status": r["final_status"],
                            "f1": float(r["f1"]) if done else None,
                            "identity": float(r["identity_pct"]) if done else None}
        return out
    for r in csv.DictReader(open(os.path.join(GRIDS, METHODS[method], "grid_results.csv"), newline="")):
        if r["depth"] not in DEPTHS:
            continue
        key = _key(r)
        if "iteration" in r and key in out and int(r["iteration"]) <= out[key]["_it"]:
            continue
        error = r.get("error", "")
        done = not error
        out[key] = {"divergence": r["divergence"], "_it": int(r.get("iteration") or 0),
                    "status": "completed" if done else ("timeout" if "timed out" in error else "no_output_or_error"),
                    "f1": float(r["f1"]) if done else None,
                    "identity": float(r["identity_pct"]) if done else None}
    if method in ("NOVOPlasty", "GetOrganelle"):
        path = os.path.join(GRIDS, METHODS[method], "normalized_successful_rerun.csv")
        norm = {_key(r): r for r in csv.DictReader(open(path, newline="")) if not r["error"]}
        for key, row in out.items():
            if row["status"] == "completed":
                n = norm[key]  # every successful job has a normalized score
                row["f1"], row["identity"] = float(n["f1"]), float(n["identity_pct"])
    return out


def _q(values: list[float]) -> tuple:
    if not values:
        return ("",) * 5
    v = sorted(values)
    if len(v) == 1:
        return (v[0], v[0], v[0], v[0], v[0])
    q1, med, q3 = st.quantiles(v, n=4, method="inclusive")
    return (v[0], q1, med, q3, v[-1])


def _cluster(ref: str, tgt: str) -> tuple:
    return tuple(sorted((ref, tgt)))


def _boot_ci(cluster_means: list[float], rng: random.Random) -> tuple:
    n = len(cluster_means)
    stats = sorted(st.mean(rng.choices(cluster_means, k=n)) for _ in range(N_BOOT))
    return stats[int(0.025 * N_BOOT)], stats[int(0.975 * N_BOOT) - 1]


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    data = {m: load(m) for m in METHODS}
    asem = data["ASEM"]
    cells = sorted({(v["divergence"], k[2]) for k, v in asem.items()},
                   key=lambda c: (["same-genus", "same-family", "same-order"].index(c[0]), int(c[1])))

    status_rows, acc_rows, paired_rows, pair_rows = [], [], [], []
    rng = random.Random(SEED)
    for div, depth in cells:
        keys = [k for k, v in asem.items() if v["divergence"] == div and k[2] == depth]
        for method, rows in data.items():
            jobs = [rows[k] for k in keys if k in rows]
            counts = defaultdict(int)
            for j in jobs:
                counts[j["status"]] += 1
            status_rows.append({"divergence": div, "depth": depth, "method": method, "n_jobs": len(keys),
                                "n_present": len(jobs), **{f"n_{s}": c for s, c in counts.items()}})
            done = [j for j in jobs if j["status"] == "completed"]
            f1s, ids = [j["f1"] for j in done], [j["identity"] for j in done]
            acc_rows.append({"divergence": div, "depth": depth, "method": method, "n_completed": len(done),
                             "f1_mean": st.mean(f1s) if f1s else "", "identity_mean": st.mean(ids) if ids else "",
                             **dict(zip(("f1_min", "f1_q1", "f1_median", "f1_q3", "f1_max"), _q(f1s))),
                             **dict(zip(("id_min", "id_q1", "id_median", "id_q3", "id_max"), _q(ids)))})
            if method == "ASEM":
                continue
            paired = [(k, rows[k]["f1"] - asem[k]["f1"], rows[k]["identity"] - asem[k]["identity"])
                      for k in keys if k in rows and rows[k]["status"] == "completed"]
            if not paired:
                paired_rows.append({"divergence": div, "depth": depth, "method": method, "n_paired": 0})
                continue
            by_cluster = defaultdict(list)
            by_ordered = defaultdict(list)
            for k, d_f1, d_id in paired:
                by_cluster[_cluster(k[0], k[1])].append((d_f1, d_id))
                by_ordered[(k[0], k[1])].append((d_f1, d_id))
            c_f1 = [st.mean(x[0] for x in v) for v in by_cluster.values()]
            c_id = [st.mean(x[1] for x in v) for v in by_cluster.values()]
            row = {"divergence": div, "depth": depth, "method": method, "n_paired": len(paired),
                   "n_clusters": len(by_cluster),
                   "dF1_mean_of_clusters": st.mean(c_f1), "dID_mean_of_clusters": st.mean(c_id),
                   **dict(zip(("dF1_job_min", "dF1_job_q1", "dF1_job_median", "dF1_job_q3", "dF1_job_max"),
                              _q([p[1] for p in paired]))),
                   **dict(zip(("dID_job_min", "dID_job_q1", "dID_job_median", "dID_job_q3", "dID_job_max"),
                              _q([p[2] for p in paired]))),
                   "share_jobs_method_better_f1": sum(p[1] > 0 for p in paired) / len(paired)}
            if len(by_cluster) >= 5:
                row["dF1_ci95_lo"], row["dF1_ci95_hi"] = _boot_ci(c_f1, rng)
                row["dID_ci95_lo"], row["dID_ci95_hi"] = _boot_ci(c_id, rng)
            paired_rows.append(row)
            for (ref, tgt), v in by_ordered.items():
                pair_rows.append({"divergence": div, "depth": depth, "method": method, "reference": ref,
                                  "target": tgt, "n_reps": len(v),
                                  "dF1_mean": st.mean(x[0] for x in v), "dF1_min": min(x[0] for x in v),
                                  "dF1_max": max(x[0] for x in v),
                                  "dID_mean": st.mean(x[1] for x in v), "dID_min": min(x[1] for x in v),
                                  "dID_max": max(x[1] for x in v)})

    def write(name: str, rows: list[dict]) -> None:
        fields: list[str] = []
        for r in rows:
            fields += [f for f in r if f not in fields]
        with open(os.path.join(OUT, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    write("status_by_cell.csv", status_rows)
    write("accuracy_by_cell.csv", acc_rows)
    write("paired_vs_asem.csv", paired_rows)
    write("pair_level.csv", pair_rows)
    print(f"wrote {len(status_rows)} status, {len(acc_rows)} accuracy, {len(paired_rows)} paired, "
          f"{len(pair_rows)} pair-level rows to {OUT}")


if __name__ == "__main__":
    main()
