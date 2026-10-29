"""Build the S2-S6 Table files (one .xlsx each: a README sheet and a data sheet).

S2  paired differences from ASEM          analysis/metrics_uncertainty/results/paired_vs_asem.csv
S3  job-level outcome of every MIA run    evaluator_n_fix/results/mia_grid_final.csv
                                          + follow-up status from mia_followup/results/mia_grid_with_followup.csv
S4  neutral-N rescoring verification      evaluator_n_fix/results/verification.csv
S5  read-length experiment, job level     read_length_all_methods/results/results.csv with the
                                          MIA/Pilon/NOVOPlasty rows replaced by the N-safe rerun
S6  empirical NA07000 reads, job level    real_reads_all_methods/results/results.csv, same replacement

S5 and S6 use the load() of the scripts that build the main-text tables, so the
files contain exactly the rows behind Tables 4 and 6. Machine names, local
output paths, and columns that are empty in every row are dropped. Every
remaining column must have a description below, or the build fails.

Run with the genome environment (needs openpyxl).
"""

from __future__ import annotations

import csv
import os
import sys

from openpyxl import Workbook
from openpyxl.cell.cell import ILLEGAL_CHARACTERS_RE
from openpyxl.styles import Font

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.join(HERE, "..", "analysis")
DROP = {"host", "theta_fasta", "assembly_fasta", "experiment", "condition"}

COMMON = {
    "method": "Method name (ASEM_recursive: optional recursive alignment, S1 Text).",
    "reference": "Starting reference sequence.",
    "target": "Target genome whose reads were simulated and against which the output was scored.",
    "divergence": "Taxonomic relationship of reference and target: same-genus, same-family, or same-order.",
    "depth": "Mean mitochondrial-target depth (X).",
    "replicate": "Technical replicate (0-2); each replicate uses its own deterministic read seed.",
    "status": "Job outcome: completed (scoreable sequence), no_output, error, or timeout.",
    "error": "Error message or tool output tail for jobs without a sequence.",
    "runtime_s": "Wall time of the job in seconds, from runs sharing the machine; not comparable across methods (see the controlled benchmark, Table 6).",
    "n_reads": "Number of reads supplied to the method.",
    "theta_len": "Length of the output sequence (bp).",
    "n_ambiguous": "Number of ambiguous bases (N or other IUPAC codes) in the output.",
    "f1": "F1 against the target (main-text definition).",
    "identity_pct": "Identity (%) against the target: correct target positions / alignment columns.",
    "recall": "Recall c/(c+u).",
    "precision": "Precision c/(c+m).",
    "c_theta": "Target positions matched by the correct base (c).",
    "w_theta": "Target positions aligned to a different or ambiguous base (m).",
    "u_t": "Target positions aligned to a gap in the output (u).",
    "theta_sha256": "SHA-256 checksum of the stored output sequence.",
    "source": "original_run, or neutral_N_rerun for MIA, Pilon, and NOVOPlasty rows rerun with retained outputs and neutral N scoring (S1 Text).",
}

TABLES = {
    "S2_Table.xlsx": (
        "S2 Table. Paired differences from ASEM on the main grid.",
        "Each row is one method, divergence tier, and depth. Differences are method minus ASEM on the same job "
        "(same reads). Clusters are unordered species pairs; the same-order tier has 13. Bootstrap intervals are "
        "given only for same-order cells with at least 5 clusters. NOVOPlasty and "
        "GetOrganelle are scored after strand and circular-origin normalization.",
        {
            "n_paired": "Jobs in which both the method and ASEM returned a sequence.",
            "n_clusters": "Unordered species-pair clusters contributing paired jobs.",
            "dF1_mean_of_clusters": "Mean over clusters of the cluster-mean F1 difference.",
            "dID_mean_of_clusters": "Mean over clusters of the cluster-mean identity difference (percentage points).",
            **{f"dF1_job_{q}": f"Job-level F1 difference: {d}." for q, d in
               [("min", "minimum"), ("q1", "first quartile"), ("median", "median"), ("q3", "third quartile"), ("max", "maximum")]},
            **{f"dID_job_{q}": f"Job-level identity difference (percentage points): {d}." for q, d in
               [("min", "minimum"), ("q1", "first quartile"), ("median", "median"), ("q3", "third quartile"), ("max", "maximum")]},
            "share_jobs_method_better_f1": "Fraction of paired jobs in which the method had the higher F1.",
            "dF1_ci95_lo": "Same-order tier only: lower bound of the 95% cluster-bootstrap interval of dF1 (10,000 resamples).",
            "dF1_ci95_hi": "Same-order tier only: upper bound of the 95% cluster-bootstrap interval of dF1.",
            "dID_ci95_lo": "Same-order tier only: lower bound of the 95% cluster-bootstrap interval of the identity difference.",
            "dID_ci95_hi": "Same-order tier only: upper bound of the 95% cluster-bootstrap interval of the identity difference.",
        }),
    "S3_Table.xlsx": (
        "S3 Table. Job-level outcome of every MIA run.",
        "All 720 main-grid MIA jobs: the original grid run (300-s limit), the follow-up of limit hits and process "
        "errors (3,600-s limit), and the final outcome after the rerun with retained outputs and neutral N scoring "
        "plus one sequential retry of failed jobs. Metrics are those reported in the main text.",
        {
            "original_status": "Outcome in the original grid run (300-s limit).",
            "original_runtime_s": "Wall time in the original grid run (s).",
            "followup_budget_s": "Time limit of the follow-up run (s); empty if the job was not followed up.",
            "followup_status": "Outcome after the 3,600-s follow-up.",
            "followup_runtime_s": "Wall time of the run that defined followup_status (s).",
            "final_status": "Final outcome after the neutral-N rerun and sequential retry.",
            "final_source": "nsafe_rerun (parallel rerun) or nsafe_retry (sequential retry of a failed job).",
        }),
    "S4_Table.xlsx": (
        "S4 Table. Rescoring of external-tool outputs with the neutral N score.",
        "One row per rerun job. Outputs with ambiguous bases were rescored with the original evaluator "
        "(N scored as a mismatch during alignment) and with the neutral N score.",
        {
            "experiment": "grid (main grid), macaca, read_length, or real_reads.",
            "condition": "Depth (grid, macaca, real_reads) or read-length condition (read_length).",
            "recorded_f1": "F1 recorded in the original run (empty if that run returned no sequence).",
            "old_evaluator_f1": "F1 of the rerun output under the original evaluator.",
            "fixed_f1": "F1 of the rerun output under the neutral N score (reported value; real_reads: phase-3 truth, primary mask).",
            "verdict": "f1_matched: the rerun output's F1 under the original evaluator equals the recorded F1, and "
                       "the neutral N score leaves it unchanged; f1_matched_rescored: equal under the original "
                       "evaluator, changed by the neutral N score; f1_not_matched: the original-evaluator F1 differs "
                       "from the recorded F1; recorded_failed: the original run had no output; rerun_failed / "
                       "both_failed: no output in the rerun / in both runs. Only F1 values are compared, because the "
                       "original runs did not retain their output sequences.",
            "target": "Target genome for simulated experiments; empty for real_reads rows, which are scored against "
                      "the curated phase-3 NA07000 sequence with the primary mask.",
        }),
    "S5_Table.xlsx": (
        "S5 Table. Read-length experiment, all methods.",
        "One representative pair per divergence level at 4X, three replicates, read lengths 50-300 bp and a "
        "variable-length condition. Cell means are in Table 4.",
        {
            "length_condition": "Read length in bp, or variable (150 bp kept with probability 0.7, otherwise trimmed to 36-149 bp).",
            "mean_read_len": "Mean length of the supplied reads (bp).",
        }),
    "S6_Table.xlsx": (
        "S6 Table. Empirical-read experiment (NA07000), all methods.",
        "Reads sampled from the mtDNA-mapped reads of NA07000 at 1-8X, three replicates, two starting references. "
        "Each output is scored against the curated phase-3 NA07000 sequence with the primary mask (reported), the "
        "strict mask, and no mask, and against rCRS. Depth means are in Table 5.",
        {
            **{f"{p}_{m}": f"{label}: {desc}" for p, label in
               [("truth_primary", "Phase-3 truth, primary mask (reported)"), ("truth_strict", "Phase-3 truth, strict mask"),
                ("truth_unmasked", "Phase-3 truth, no mask"), ("rcrs", "rCRS")]
               for m, desc in [("identity_pct", "identity (%)."), ("recall", "recall."), ("precision", "precision."),
                               ("f1", "F1."), ("c_theta", "correct target positions (c)."),
                               ("w_theta", "incorrect target positions (m)."), ("u_t", "unrecovered target positions (u).")]},
        }),
}


def read(*parts: str) -> list[dict]:
    with open(os.path.join(ANALYSIS, *parts), newline="") as fh:
        return list(csv.DictReader(fh))


def load_with_overlay(module_dir: str, module: str, rerun_file: str) -> list[dict]:
    sys.path.insert(0, os.path.join(ANALYSIS, module_dir))
    rows = __import__(module).load()
    fixed = {r["method"] for r in read("evaluator_n_fix", "results", rerun_file)}
    for r in rows:
        r["source"] = "neutral_N_rerun" if r["method"] in fixed else "original_run"
    return rows


def write(name: str, rows: list[dict], keep_experiment: bool = False) -> None:
    title, note, specific = TABLES[name]
    drop = DROP - ({"experiment", "condition"} if keep_experiment else set())
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields and k not in drop]
    fields = [k for k in fields if any(r.get(k) not in ("", None) for r in rows)]
    desc = {**COMMON, **specific}
    missing = [k for k in fields if k not in desc]
    if missing:
        raise SystemExit(f"{name}: no description for columns {missing}")
    wb = Workbook()
    ws = wb.active
    ws.title = "README"
    ws.append([title])
    ws["A1"].font = Font(bold=True)
    ws.append([note])
    ws.append([])
    ws.append(["Column", "Description"])
    for cell in ws[4]:
        cell.font = Font(bold=True)
    for k in fields:
        ws.append([k, desc[k]])
    ws.column_dimensions["A"].width = 30
    ws.column_dimensions["B"].width = 110
    data = wb.create_sheet("data")
    data.append(fields)
    for cell in data[1]:
        cell.font = Font(bold=True)
    for r in rows:
        data.append([_value(r.get(k, "")) for k in fields])
    data.freeze_panes = "A2"
    wb.save(os.path.join(HERE, name))
    print(f"{name}: {len(rows)} rows, {len(fields)} columns")


def _value(v):
    """Numbers as numbers, everything else as text."""
    if isinstance(v, (int, float)) or v in ("", None):
        return v
    try:
        return int(v) if v.lstrip("-").isdigit() else float(v)
    except ValueError:
        # Tool logs carry control characters that Excel rejects; keep a short, clean tail.
        v = ILLEGAL_CHARACTERS_RE.sub(" ", v)
        return v if len(v) <= 300 else v[:300] + " ..."


def main() -> None:
    write("S2_Table.xlsx", read("metrics_uncertainty", "results", "paired_vs_asem.csv"))

    key = lambda r: (r["reference"], r["target"], r["depth"], r["replicate"])
    followup = {key(r): r for r in read("mia_followup", "results", "mia_grid_with_followup.csv")}
    s3 = []
    for r in read("evaluator_n_fix", "results", "mia_grid_final.csv"):
        f = followup[key(r)]
        row = {k: r[k] for k in ("reference", "target", "divergence", "depth", "replicate",
                                 "original_status", "original_runtime_s", "followup_budget_s")}
        row["followup_status"] = f["final_status"]
        row["followup_runtime_s"] = f["final_runtime_s"]
        row["final_status"] = r["final_status"]
        row["final_source"] = r["source"]
        row.update({k: r[k] for k in ("theta_len", "n_ambiguous", "f1", "identity_pct",
                                      "recall", "precision", "theta_sha256")})
        s3.append(row)
    write("S3_Table.xlsx", s3)

    write("S4_Table.xlsx", read("evaluator_n_fix", "results", "verification.csv"), keep_experiment=True)
    write("S5_Table.xlsx", load_with_overlay("read_length_all_methods", "make_table_read_length", "read_length.csv"))
    write("S6_Table.xlsx", load_with_overlay("real_reads_all_methods", "make_table_real_reads", "real_reads.csv"))


if __name__ == "__main__":
    main()
