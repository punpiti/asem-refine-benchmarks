"""Check the reruns against the originally recorded values.

For every rerun row:
  * N-free output: nsafe falls back to the original aligner call, so its
    metrics under both rules are the same.
  * Output containing N: the sequence is re-scored with the ORIGINAL
    (unpatched) evaluator.
  The rerun F1 under the original evaluator is then compared with the
  recorded F1. Only F1 values are compared: the original runs did not retain
  their output sequences, so equal F1 is consistent with, but does not prove,
  an identical sequence.

Verdicts: f1_matched (original-evaluator F1 equals the recorded F1 and the
neutral rule leaves it unchanged), f1_matched_rescored (equal under the
original evaluator, changed by the neutral rule), f1_not_matched,
recorded_failed, rerun_failed, both_failed, no_record.

Writes results/verification.csv and prints a summary.
"""

from __future__ import annotations

import csv
import gzip
import os
import sys
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(HERE)
BENCH_SCRIPTS = os.environ.get("BENCH_SCRIPTS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "scripts"))
sys.path.insert(0, BENCH_SCRIPTS)
sys.path.insert(0, os.path.join(ANALYSIS, "real_reads_all_methods"))
sys.path.insert(0, os.path.join(ANALYSIS, "read_length_all_methods"))

import common  # noqa: E402  (unpatched)
from common import read_fasta  # noqa: E402

RES = os.path.join(HERE, "results")
B = os.path.join(BENCH_SCRIPTS, "results")
TOL = 1e-9


def _load(path, key):
    with open(path, newline="") as fh:
        return {key(r): r for r in csv.DictReader(fh)}


def recorded() -> dict:
    """(experiment, method, ref, target, condition, replicate) -> recorded F1 string."""
    rec = {}
    g = lambda r: (r["reference"], r["target"], str(int(r["depth"])), str(int(r["replicate"])))
    mia = _load(os.path.join(ANALYSIS, "mia_followup", "results", "mia_grid_with_followup.csv"), g)
    for k, r in mia.items():
        rec[("grid", "MIA", *k)] = r["f1"] if r["final_status"] == "completed" else None
    for k, r in _load(os.path.join(B, "ext_pilon", "grid_results.csv"), g).items():
        rec[("grid", "Pilon", *k)] = r["f1"] if not r["error"] else None
    for m in ("MIA", "Pilon", "NOVOPlasty", "GetOrganelle"):
        for k, r in _load(os.path.join(B, "same_genus_macaca", f"ext_{m.lower()}", "grid_results.csv"), g).items():
            rec[("macaca", m, *k)] = r["f1"] if not r["error"] else None
    rl = lambda r: (r["method"], r["reference"], r["target"], r["length_condition"], str(int(r["replicate"])))
    for k, r in _load(os.path.join(ANALYSIS, "read_length_all_methods", "results", "results.csv"), rl).items():
        rec[("read_length", *k)] = r["f1"] if r["status"] == "completed" else None
    rr = lambda r: (r["method"], r["reference"], "", str(int(r["depth"])), str(int(r["replicate"])))
    for k, r in _load(os.path.join(ANALYSIS, "real_reads_all_methods", "results", "results.csv"), rr).items():
        rec[("real_reads", *k)] = r["truth_primary_f1"] if r["status"] == "completed" else None
    return rec


def _key(r: dict) -> tuple:
    cond = r["condition"] or r["length_condition"] or r["depth"]
    if r["experiment"] in ("grid", "macaca", "real_reads"):
        cond = str(int(float(cond)))
    return (r["experiment"], r["method"], r["reference"], r.get("target", ""), cond, str(int(r["replicate"])))


def _old_f1(row: dict) -> float:
    with gzip.open(os.path.join(HERE, row["theta_fasta"]), "rt") as fh:
        theta = fh.read().split("\n")[1]
    if row["experiment"] == "real_reads":
        import run_real_reads_all_methods as mod  # unpatched module-level aligner
        truth = read_fasta(mod.TRUTH_FASTA)[1].replace("-", "")
        mask = {i for p in mod.PRIMARY_MASK_RCRS if (i := mod.rcrs_to_truth_index(p)) is not None} \
            if not hasattr(mod, "PRIMARY_MASK") else mod.PRIMARY_MASK
        return mod.masked_metrics(theta, truth, mask, True)["f1"]
    # Pool workers are reused across rows, so set the pair set for every row.
    os.environ["ASEM_PAIR_SET"] = "macaca_same_genus" if row["experiment"] == "macaca" else "phase1"
    import importlib
    import phase1_species
    importlib.reload(phase1_species)
    target = read_fasta(phase1_species.fasta_path(row["target"]))[1]
    return common.evaluate_theta_vs_target(theta, target, normalize_strand=True,
                                           normalize_circular_origin=True)["f1"]


def main() -> None:
    rec = recorded()
    rows = []
    for exp in ("grid", "macaca", "read_length", "real_reads"):
        path = os.path.join(RES, f"{exp}.csv")
        if os.path.exists(path):
            rows += list(csv.DictReader(open(path, newline="")))
    # A failed grid job that completed in the sequential retry is verified
    # through its retry output (see merge_mia_grid.py).
    retry_path = os.path.join(RES, "grid_crash_retry.csv")
    if os.path.exists(retry_path):
        done = {_key(r): r for r in csv.DictReader(open(retry_path, newline="")) if r["status"] == "completed"}
        rows = [done.get(_key(r), r) if r["experiment"] == "grid" else r for r in rows]
    need_old = [r for r in rows if r["status"] == "completed" and int(r["n_ambiguous"] or 0) > 0]
    # Re-scoring is slow (quadratic-memory alignment), so results are cached
    # by output file and a failed run resumes where it stopped.
    cache_path = os.path.join(RES, "verification_old_f1_cache.csv")
    cache = {}
    if os.path.exists(cache_path):
        cache = {r["theta_fasta"]: float(r["old_f1"]) for r in csv.DictReader(open(cache_path, newline=""))}
    todo = [r for r in need_old if r["theta_fasta"] not in cache]
    is_new = not os.path.exists(cache_path)
    with open(cache_path, "a", newline="") as fh, \
            ProcessPoolExecutor(max_workers=int(os.environ.get("VERIFY_PARALLEL", "8"))) as pool:
        w = csv.writer(fh)
        if is_new:
            w.writerow(["theta_fasta", "old_f1"])
        for r, f1 in zip(todo, pool.map(_old_f1, todo)):
            cache[r["theta_fasta"]] = f1
            w.writerow([r["theta_fasta"], f1])
            fh.flush()
    old = {id(r): cache[r["theta_fasta"]] for r in need_old}
    out = []
    for r in rows:
        k = _key(r)
        recorded_f1 = rec.get(k, "missing")
        new_f1 = r.get("f1") or r.get("truth_primary_f1") or ""
        if r["status"] != "completed":
            verdict = "both_failed" if recorded_f1 is None else "rerun_failed"
            old_f1 = ""
        else:
            old_f1 = old.get(id(r), float(new_f1))
            if recorded_f1 in (None, "missing"):
                verdict = "recorded_failed" if recorded_f1 is None else "no_record"
            elif abs(float(recorded_f1) - old_f1) <= TOL:
                verdict = "f1_matched" if abs(old_f1 - float(new_f1)) <= TOL else "f1_matched_rescored"
            else:
                verdict = "f1_not_matched"
        out.append({"experiment": k[0], "method": k[1], "reference": k[2], "target": k[3], "condition": k[4],
                    "replicate": k[5], "n_ambiguous": r.get("n_ambiguous", ""), "recorded_f1": recorded_f1,
                    "old_evaluator_f1": old_f1, "fixed_f1": new_f1, "verdict": verdict})
    with open(os.path.join(RES, "verification.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    from collections import Counter
    for (exp, m, v), n in sorted(Counter((o["experiment"], o["method"], o["verdict"]) for o in out).items()):
        print(f"{exp:12s}{m:13s}{v:26s}{n}")


if __name__ == "__main__":
    main()
