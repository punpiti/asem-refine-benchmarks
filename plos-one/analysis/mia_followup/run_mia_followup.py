"""Follow-up of censored and crashed MIA jobs from the v1.3.0 full grid.

The released grid (asem-refine-benchmarks v1.3.0,
scripts/results/ext_mia/grid_results.csv) stopped every MIA run at 300 s and
recorded 69 timeouts plus 18 process errors (return codes -6 and -11). This
runner re-executes exactly those job keys with the original reference/target
sequences, read simulator, stable read seed, MIA flags and evaluator, under a
larger, explicit time budget. It never touches the released CSV.

Each attempt records the original status, the new status, the budget, wall
time, peak RSS, return code, retained stdout/stderr logs and (on success) the
consensus FASTA with its SHA-256. Results are appended to a separate CSV and
the runner is resumable per (job key, budget).

Usage (from the project root):

    micromamba run -n genome python plos-one/analysis/mia_followup/run_mia_followup.py \
        --select pilot --timeout 3600 --parallel 6
    ... --select timeouts --timeout 3600 --parallel 6
    ... --select crashes --timeout 3600 --parallel 6
"""

from __future__ import annotations

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse  # noqa: E402
import csv  # noqa: E402
import fcntl  # noqa: E402
import glob  # noqa: E402
import gzip  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import platform  # noqa: E402
import shutil  # noqa: E402
import signal  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from datetime import datetime, timezone  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BENCH_REPO = PROJECT_ROOT
BENCH_SCRIPTS = os.path.join(BENCH_REPO, "scripts")
sys.path.insert(0, BENCH_SCRIPTS)

from common import evaluate_theta_vs_target, read_fasta, simulate_shotgun_reads, stable_seed  # noqa: E402
from ext_mia import MA_BIN, MIA_BIN, _parse_single_fasta, _write_fasta  # noqa: E402
from phase1_species import fasta_path  # noqa: E402

ORIGINAL_CSV = os.path.join(BENCH_SCRIPTS, "results", "ext_mia", "grid_results.csv")
RESULTS_DIR = os.path.join(HERE, "results")
FOLLOWUP_CSV = os.path.join(RESULTS_DIR, "followup_results.csv")
LOG_DIR = os.path.join(RESULTS_DIR, "logs")
ASSEMBLY_DIR = os.path.join(RESULTS_DIR, "assemblies")
MANIFEST_DIR = os.path.join(RESULTS_DIR, "manifests")
READ_LENGTH = 150

FIELDNAMES = [
    "reference", "target", "divergence", "depth", "replicate", "n_reads",
    "original_status", "original_error", "original_runtime_s",
    "timeout_budget_s", "status", "returncode", "wall_s", "peak_rss_kb",
    "output_present", "theta_len", "n_ambiguous",
    "identity_pct", "recall", "precision", "f1", "c_theta", "w_theta", "u_t",
    "stdout_log", "stderr_log", "assembly_fasta", "assembly_sha256",
    "error", "started_utc", "host",
]

_SEQ_CACHE: dict[str, str] = {}


def _seq(name: str) -> str:
    if name not in _SEQ_CACHE:
        _SEQ_CACHE[name] = read_fasta(fasta_path(name))[1]
    return _SEQ_CACHE[name]


def _sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _original_status(row: dict) -> str:
    error = row["error"]
    if not error:
        return "completed"
    if error.startswith("timed out"):
        return "timeout"
    if error.startswith("MIA exited"):
        return "process_error"
    return "other_error"


def load_original() -> list[dict]:
    with open(ORIGINAL_CSV, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        row["original_status"] = _original_status(row)
    return rows


def select_jobs(rows: list[dict], selection: str) -> list[dict]:
    timeouts = [r for r in rows if r["original_status"] == "timeout"]
    crashes = [r for r in rows if r["original_status"] == "process_error"]
    if selection == "timeouts":
        return timeouts
    if selection == "crashes":
        return crashes
    if selection == "pilot":
        # One timed-out job per (divergence, depth) pattern: the first by key order.
        chosen: dict[tuple[str, str], dict] = {}
        for row in sorted(timeouts, key=_key):
            chosen.setdefault((row["divergence"], row["depth"]), row)
        return sorted(chosen.values(), key=lambda r: (r["divergence"], int(r["depth"])))
    raise ValueError(selection)


def _key(row: dict) -> tuple[str, str, int, int]:
    return (row["reference"], row["target"], int(row["depth"]), int(row["replicate"]))


def _key_slug(row: dict) -> str:
    return f"{row['reference']}__{row['target']}__{int(row['depth'])}X__r{int(row['replicate'])}"


def _gzip_text(text: str, path: str) -> None:
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def run_one(row: dict, timeout: int) -> dict:
    reference = _seq(row["reference"])
    target = _seq(row["target"])
    depth, replicate = int(row["depth"]), int(row["replicate"])
    reads = simulate_shotgun_reads(target, float(depth), READ_LENGTH, stable_seed(row["target"], depth, replicate))
    slug = f"{_key_slug(row)}__t{timeout}"
    out = {
        "reference": row["reference"], "target": row["target"],
        "divergence": row["divergence"], "depth": depth, "replicate": replicate,
        "n_reads": len(reads),
        "original_status": row["original_status"],
        "original_error": row["error"].splitlines()[0][:80] if row["error"] else "",
        "original_runtime_s": row["runtime_s"],
        "timeout_budget_s": timeout,
        "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "host": platform.node(),
        "output_present": False,
        "theta_len": 0, "n_ambiguous": 0,
        "error": "",
    }
    if len(reads) != int(row["n_reads"]):
        # Read regeneration must reproduce the original input exactly.
        return {**out, "status": "input_mismatch", "returncode": "", "wall_s": 0,
                "error": f"n_reads {len(reads)} != original {row['n_reads']}"}

    workdir = tempfile.mkdtemp(prefix="mia_followup_")
    stdout_log = os.path.join(LOG_DIR, f"{slug}.stdout.gz")
    stderr_log = os.path.join(LOG_DIR, f"{slug}.stderr.gz")
    try:
        reference_path = os.path.join(workdir, "reference.fasta")
        reads_path = os.path.join(workdir, "reads.fasta")
        output_root = os.path.join(workdir, "assembly")
        time_path = os.path.join(workdir, "time.txt")
        _write_fasta(reference_path, [("reference", reference)])
        _write_fasta(reads_path, [(f"read{i}", read) for i, read in enumerate(reads)])
        cmd = [
            "/usr/bin/time", "-f", "%e %M", "-o", time_path,
            MIA_BIN, "-r", reference_path, "-f", reads_path, "-m", output_root,
            "-c", "-D", "-F",
        ]
        started = time.time()
        proc = subprocess.Popen(
            cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, start_new_session=True,
        )
        try:
            stdout, stderr = proc.communicate(timeout=timeout)
            timed_out = False
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            stdout, stderr = proc.communicate()
            timed_out = True
        wall = time.time() - started
        _gzip_text(stdout, stdout_log)
        _gzip_text(stderr, stderr_log)
        out.update(stdout_log=os.path.relpath(stdout_log, HERE),
                   stderr_log=os.path.relpath(stderr_log, HERE), wall_s=round(wall, 2))
        peak_rss = ""
        if os.path.exists(time_path):
            # /usr/bin/time prints a "Command terminated by signal" line first on abnormal exit.
            fields = open(time_path).read().strip().splitlines()[-1].split()
            if len(fields) == 2 and fields[1].isdigit():
                peak_rss = int(fields[1])
        out["peak_rss_kb"] = peak_rss

        # /usr/bin/time reports a child killed by signal N as exit 128+N.
        rc = proc.returncode
        if not timed_out and rc is not None and rc > 128:
            rc = -(rc - 128)
        out["returncode"] = "" if timed_out else rc
        if timed_out:
            return {**out, "status": "timeout", "error": f"timed out after {timeout}s"}
        if rc != 0:
            return {**out, "status": "process_error", "error": f"MIA exited {rc}; stderr tail: {stderr[-300:]}"}

        maln = glob.glob(f"{output_root}.*")
        if not maln:
            return {**out, "status": "no_output", "error": "MIA produced no final MALN assembly"}
        final_maln = max(maln, key=lambda path: int(path.rsplit(".", 1)[1]))
        converted = subprocess.run(
            [MA_BIN, "-M", final_maln, "-f", "5", "-I", "MIA_consensus"],
            cwd=workdir, capture_output=True, text=True, timeout=600,
        )
        if converted.returncode != 0:
            return {**out, "status": "conversion_error", "error": f"ma exited {converted.returncode}"}
        try:
            theta = _parse_single_fasta(converted.stdout)
        except Exception as exc:  # empty consensus
            return {**out, "status": "empty_output", "error": str(exc)}
        assembly = os.path.join(ASSEMBLY_DIR, f"{slug}.fasta")
        _write_fasta(assembly, [(f"MIA_consensus {slug}", theta)])
        metrics = evaluate_theta_vs_target(theta, target, normalize_strand=True, normalize_circular_origin=True)
        return {
            **out, **metrics, "status": "completed", "output_present": True,
            "theta_len": len(theta), "n_ambiguous": sum(s not in "ACGT" for s in theta),
            "assembly_fasta": os.path.relpath(assembly, HERE), "assembly_sha256": _sha256(assembly),
        }
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def _done_keys() -> set[tuple]:
    if not os.path.exists(FOLLOWUP_CSV):
        return set()
    with open(FOLLOWUP_CSV, newline="") as fh:
        return {(*_key(r), int(r["timeout_budget_s"])) for r in csv.DictReader(fh)}


def write_manifest(args: argparse.Namespace, jobs: list[dict]) -> str:
    def _git(*cmd: str) -> str:
        try:
            return subprocess.run(["git", "-C", BENCH_REPO, *cmd], capture_output=True, text=True).stdout.strip()
        except OSError:
            return ""

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "host": platform.node(),
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "python": sys.version,
        "selection": args.select,
        "timeout_s": args.timeout,
        "parallel": args.parallel,
        "n_jobs": len(jobs),
        "jobs": [list(_key(r)) for r in jobs],
        "mia_bin": MIA_BIN, "mia_sha256": _sha256(MIA_BIN),
        "ma_bin": MA_BIN, "ma_sha256": _sha256(MA_BIN),
        "mia_package": "bioconda mapping-iterative-assembler 1.0 h503566f_7",
        "mia_flags": ["-c", "-D", "-F"],
        "read_length": READ_LENGTH,
        "benchmark_repo_commit": _git("rev-parse", "HEAD"),
        "benchmark_repo_describe": _git("describe", "--tags", "--always"),
        "benchmark_repo_dirty_scripts": _git("status", "--porcelain", "--", "scripts/common.py",
                                             "scripts/ext_mia.py", "scripts/phase1_species.py"),
        "original_csv": os.path.relpath(ORIGINAL_CSV, PROJECT_ROOT),
        "original_csv_sha256": _sha256(ORIGINAL_CSV),
    }
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    path = os.path.join(MANIFEST_DIR, f"manifest_{args.select}_t{args.timeout}_{stamp}.json")
    with open(path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--select", choices=["pilot", "timeouts", "crashes"], required=True)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--parallel", type=int, default=4)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    for path in (RESULTS_DIR, LOG_DIR, ASSEMBLY_DIR, MANIFEST_DIR):
        os.makedirs(path, exist_ok=True)
    lock = open(os.path.join(RESULTS_DIR, ".followup.lock"), "w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print("another follow-up runner holds the lock; exiting", flush=True)
        return

    selected = select_jobs(load_original(), args.select)
    done = _done_keys()
    jobs = [r for r in selected if (*_key(r), args.timeout) not in done]
    print(f"=== MIA follow-up [{args.select}] budget={args.timeout}s: "
          f"{len(selected)} selected, {len(selected) - len(jobs)} already done, "
          f"{len(jobs)} to run, {args.parallel} parallel ===", flush=True)
    for row in jobs:
        print("  ", _key_slug(row), row["divergence"], row["original_status"], flush=True)
    if args.dry_run or not jobs:
        return
    print("manifest:", write_manifest(args, jobs), flush=True)

    is_new = not os.path.exists(FOLLOWUP_CSV)
    with open(FOLLOWUP_CSV, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        if is_new:
            writer.writeheader()
            fh.flush()
        with ProcessPoolExecutor(max_workers=args.parallel) as pool:
            futures = {pool.submit(run_one, row, args.timeout): row for row in jobs}
            for n, future in enumerate(as_completed(futures), start=1):
                result = future.result()
                writer.writerow({k: result.get(k, "") for k in FIELDNAMES})
                fh.flush()
                print(f"[{n}/{len(jobs)}] {result['reference']} -> {result['target']} "
                      f"{result['depth']}X r{result['replicate']}: {result['status']} "
                      f"wall={result.get('wall_s', '')}s rss={result.get('peak_rss_kb', '')}KB "
                      f"F1={result.get('f1', '')}", flush=True)


if __name__ == "__main__":
    main()
