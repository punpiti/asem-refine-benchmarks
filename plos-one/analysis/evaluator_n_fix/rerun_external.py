"""Rerun external tools whose output sequences were not retained, keep every
sequence, and score it with the N-safe evaluator (nsafe.py).

Experiments (job keys and read generation exactly as in the original runs):
  grid        MIA and Pilon, 720 phase-1 keys each (reads: stable_seed(target,
              depth, replicate), 150 bp); MIA limit 3600 s as in the follow-up.
  macaca      MIA, Pilon, NOVOPlasty, GetOrganelle on the Macaca pair
              (run with ASEM_PAIR_SET=macaca_same_genus).
  read_length MIA, Pilon, NOVOPlasty through run_read_length_all_methods.run_one.
  real_reads  MIA, Pilon, NOVOPlasty through run_real_reads_all_methods.run_one.

Each output sequence is written to results/thetas/<experiment>/<method>/ and
its SHA-256 recorded; verify_unchanged.py later re-scores the same sequences
with the original evaluator and compares with the recorded values.

    ASEM_PAIR_SET=... python -I rerun_external.py --experiment grid --parallel 12
"""

from __future__ import annotations

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse  # noqa: E402
import csv  # noqa: E402
import fcntl  # noqa: E402
import gzip  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import platform  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from datetime import datetime, timezone  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ANALYSIS, "read_length_all_methods"))
sys.path.insert(0, os.path.join(ANALYSIS, "real_reads_all_methods"))

import nsafe  # noqa: E402
import common  # noqa: E402
from common import read_fasta, simulate_shotgun_reads, stable_seed  # noqa: E402
from phase1_species import fasta_path  # noqa: E402
import ext_getorganelle  # noqa: E402
import ext_novoplasty  # noqa: E402
from ext_mia import run_mia  # noqa: E402
from ext_pilon import run_pilon  # noqa: E402

nsafe.install()

BENCH_RESULTS = os.path.join(nsafe.BENCH_SCRIPTS, "results")
# Workers start via forkserver and re-import this module, so the smoke-test
# redirection travels in the environment rather than in module globals.
RESULTS_DIR = os.environ.get("NFIX_RESULTS_DIR", os.path.join(HERE, "results"))
THETA_DIR = os.path.join(RESULTS_DIR, "thetas")
TIMEOUT_S = 3600

GRID_SOURCES = {  # method -> original CSV defining the job keys
    "grid": {"MIA": "ext_mia/grid_results.csv", "Pilon": "ext_pilon/grid_results.csv"},
    "macaca": {m: f"same_genus_macaca/ext_{m.lower()}/grid_results.csv"
               for m in ("MIA", "Pilon", "NOVOPlasty", "GetOrganelle")},
}

_SEQ: dict[str, str] = {}
_CAPTURED: dict[str, str] = {}


def _seq(name: str) -> str:
    if name not in _SEQ:
        _SEQ[name] = read_fasta(fasta_path(name))[1]
    return _SEQ[name]


def save_theta(experiment: str, method: str, slug: str, theta: str) -> tuple[str, str]:
    path = os.path.join(THETA_DIR, experiment, method, f"{slug}.fa.gz")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = f">{method} {slug}\n{theta}\n".encode()
    with gzip.GzipFile(path, "wb", mtime=0) as fh:
        fh.write(data)
    return os.path.relpath(path, HERE), hashlib.sha256(theta.encode()).hexdigest()


def _call(method: str, reference: str, reads: list[str], seed: int, project: str) -> str:
    if method == "MIA":
        return run_mia(reference, reads, timeout=TIMEOUT_S)
    if method == "Pilon":
        return run_pilon(reference, reads, timeout=TIMEOUT_S)
    if method == "NOVOPlasty":
        return ext_novoplasty.run_novoplasty(reference, reads, project=project, timeout=TIMEOUT_S)
    if method == "GetOrganelle":
        return ext_getorganelle.run_getorganelle(reference, reads, quality_seed=seed, timeout=TIMEOUT_S)
    raise ValueError(method)


def grid_job(job: tuple) -> dict:
    experiment, method, ref_name, target_name, depth, replicate, n_reads = job
    target = _seq(target_name)
    seed = stable_seed(target_name, depth, replicate)
    reads = simulate_shotgun_reads(target, float(depth), 150, seed)
    slug = f"{ref_name}__{target_name}__{depth}X__r{replicate}"
    row = {"experiment": experiment, "method": method, "reference": ref_name, "target": target_name,
           "condition": depth, "replicate": replicate, "n_reads": len(reads), "host": platform.node()}
    if n_reads and len(reads) != int(n_reads):
        return {**row, "status": "input_mismatch", "error": f"n_reads {len(reads)} != {n_reads}"}
    started = time.time()
    try:
        theta = _call(method, _seq(ref_name), reads, seed, f"nfix_{slug}")
    except Exception as exc:
        msg = str(exc)
        status = "timeout" if "timed out" in msg.lower() else "no_output_or_error"
        return {**row, "status": status, "runtime_s": round(time.time() - started, 2), "error": msg[:300]}
    runtime = time.time() - started
    path, sha = save_theta(experiment, method, slug, theta)
    scores = common.evaluate_theta_vs_target(theta, target, normalize_strand=True, normalize_circular_origin=True)
    return {**row, "status": "completed", "runtime_s": round(runtime, 2), "theta_len": len(theta),
            "n_ambiguous": sum(c not in "ACGT" for c in theta), "theta_fasta": path, "theta_sha256": sha,
            **{k: scores[k] for k in ("f1", "identity_pct", "recall", "precision", "c_theta", "w_theta", "u_t")}}


def _wrap_capture(methods_dict: dict) -> None:
    """Make each method callable stash its output so run_one's result can be paired with the sequence."""
    for name, (fn, external) in list(methods_dict.items()):
        def wrapped(*a, _fn=fn, _name=name):
            theta = _fn(*a)
            _CAPTURED["theta"] = theta
            return theta
        methods_dict[name] = (wrapped, external)


def runner_job(job: tuple) -> dict:
    experiment, inner = job[0], job[1:]
    _CAPTURED.clear()
    if experiment == "read_length":
        import run_read_length_all_methods as mod
        result = mod.run_one(inner)
        cond = str(inner[3])
        slug = f"{inner[1]}__{inner[2]}__{cond}__r{inner[4]}"
        method = inner[0]
    else:
        import run_real_reads_all_methods as mod
        result = mod.run_one(inner)
        method, ref_name, depth, replicate = inner
        slug = f"{ref_name}__{depth}X__r{replicate}"
    out = {"experiment": experiment, **result}
    if "theta" in _CAPTURED:
        path, sha = save_theta(experiment, method, slug, _CAPTURED["theta"])
        out.update(theta_fasta=path, theta_sha256=sha,
                   n_ambiguous=sum(c not in "ACGT" for c in _CAPTURED["theta"]))
    return out


def _init_runner(experiment: str) -> None:
    if experiment == "read_length":
        import run_read_length_all_methods as mod
        _wrap_capture(mod.METHODS)
    else:
        import run_real_reads_all_methods as mod
        nsafe.install(mod)
        _wrap_capture(mod.METHODS)  # same dict object as read-length METHODS
        mod._init(mod.extract_reads(), read_fasta(mod.TRUTH_FASTA)[1].replace("-", ""),
                  read_fasta(mod.RCRS_FASTA)[1])


def jobs_for(experiment: str, methods: list[str]) -> list[tuple]:
    if experiment in GRID_SOURCES:
        jobs = []
        for method in methods:
            with open(os.path.join(BENCH_RESULTS, GRID_SOURCES[experiment][method]), newline="") as fh:
                for r in csv.DictReader(fh):
                    jobs.append((experiment, method, r["reference"], r["target"], int(r["depth"]),
                                 int(r["replicate"]), r.get("n_reads", "")))
        return jobs
    if experiment == "read_length":
        import run_read_length_all_methods as mod
        return [("read_length", m, ref, tgt, lc, rep) for m in methods for ref, tgt in mod.PAIRS
                for lc in mod.LENGTH_CONDITIONS for rep in range(mod.N_REPLICATES)]
    import run_real_reads_all_methods as mod
    return [("real_reads", m, ref, d, rep) for ref in mod.REFERENCES for d in mod.DEPTHS
            for rep in range(mod.N_REPLICATES) for m in methods]


FIELDS = ["experiment", "method", "reference", "target", "divergence", "condition", "length_condition",
          "depth", "replicate", "n_reads", "status", "error", "runtime_s", "theta_len", "n_ambiguous",
          "f1", "identity_pct", "recall", "precision", "c_theta", "w_theta", "u_t",
          "theta_fasta", "theta_sha256", "host"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", choices=["grid", "macaca", "read_length", "real_reads"], required=True)
    parser.add_argument("--methods", default="")
    parser.add_argument("--parallel", type=int, default=8)
    parser.add_argument("--limit", type=int, default=0, help="smoke test: first N jobs, results under results/smoke/")
    args = parser.parse_args()
    global RESULTS_DIR, THETA_DIR
    if args.limit:
        RESULTS_DIR = os.path.join(RESULTS_DIR, "smoke")
        THETA_DIR = os.path.join(RESULTS_DIR, "thetas")
        os.environ["NFIX_RESULTS_DIR"] = RESULTS_DIR
    default = {"grid": "MIA,Pilon", "macaca": "MIA,Pilon,NOVOPlasty,GetOrganelle",
               "read_length": "MIA,Pilon,NOVOPlasty", "real_reads": "MIA,Pilon,NOVOPlasty"}
    methods = (args.methods or default[args.experiment]).split(",")
    os.makedirs(RESULTS_DIR, exist_ok=True)
    lock = open(os.path.join(RESULTS_DIR, f".{args.experiment}.lock"), "w")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

    out_csv = os.path.join(RESULTS_DIR, f"{args.experiment}.csv")
    done: set[tuple] = set()
    if os.path.exists(out_csv):
        if args.experiment not in GRID_SOURCES:
            sys.exit(f"{out_csv} exists; resume is supported only for grid experiments")
        with open(out_csv, newline="") as fh:
            done = {(r["method"], r["reference"], r["target"], int(r["condition"]), int(r["replicate"]))
                    for r in csv.DictReader(fh)}
    jobs = [j for j in jobs_for(args.experiment, methods) if tuple(j[1:6]) not in done]
    if args.limit:
        jobs = [j for i, j in enumerate(jobs) if i % max(1, len(jobs) // args.limit) == 0][:args.limit]
    print(f"=== {args.experiment}: {len(jobs)} jobs, {args.parallel} parallel ===", flush=True)
    with open(os.path.join(RESULTS_DIR, f"manifest_{args.experiment}_{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}.json"), "w") as fh:
        json.dump({"experiment": args.experiment, "methods": methods, "n_jobs": len(jobs),
                   "pair_set": os.environ.get("ASEM_PAIR_SET", "phase1"), "bench_scripts": nsafe.BENCH_SCRIPTS,
                   "host": platform.node(), "python": sys.version}, fh, indent=2)

    is_runner = args.experiment in ("read_length", "real_reads")
    fields = list(FIELDS)
    if args.experiment == "real_reads":
        import run_real_reads_all_methods as mod
        fields += [f for f in mod.run_fields() if f not in fields]
    is_new = not os.path.exists(out_csv)
    with open(out_csv, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        if is_new:
            writer.writeheader()
        pool_kwargs = {"initializer": _init_runner, "initargs": (args.experiment,)} if is_runner else {}
        with ProcessPoolExecutor(max_workers=args.parallel, **pool_kwargs) as pool:
            futures = [pool.submit(runner_job if is_runner else grid_job, j) for j in jobs]
            for n, fut in enumerate(as_completed(futures), start=1):
                row = fut.result()
                writer.writerow({k: row.get(k, "") for k in fields})
                fh.flush()
                print(f"[{n}/{len(jobs)}] {row.get('method')} {row.get('reference')} {row.get('target', '')} "
                      f"{row.get('condition', row.get('length_condition', row.get('depth')))} r{row.get('replicate')} "
                      f"{row.get('status')} f1={row.get('f1', '')}", flush=True)


if __name__ == "__main__":
    main()
