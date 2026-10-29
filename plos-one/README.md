# plos-one/: analyses added in v1.4.0

These scripts and result files produce every number, table, and data figure of
the revised ASEM/ASEM-Hybrid manuscript that is not already in
`scripts/results/`. They build on the v1.3.0 benchmark code in `scripts/`,
which is unchanged in this release.

## Scoring of ambiguous bases (neutral N)

The evaluator in `scripts/common.py` aligns each output to the target with
scikit-bio's `pair_align_nucl`, which scores N as a mismatch (-3). An output in
which a large share of bases are N (MIA marks weakly supported columns with N)
could then align best by leaving the target unaligned and was scored F1 = 0
although thousands of its bases were correct. `analysis/evaluator_n_fix/nsafe.py`
scores N (and other ambiguity codes, converted to N) as 0 during the alignment
step only; N still counts as an incorrect base in the metrics. For sequences
without ambiguous bases the original aligner call is used unchanged. All values
reported in the manuscript use this rule. `scripts/common.py` itself is not
modified, so the v1.3.0 results remain reproducible as released.

## Layout

| Folder | What it contains |
|---|---|
| `analysis/evaluator_n_fix/` | Neutral-N evaluator, rerun of every main-grid MIA job and of the external-tool jobs of the *Macaca*, read-length, and empirical-read experiments with output sequences retained (`results/thetas/`, with SHA-256 checksums), sequential retry of failed MIA jobs, verification against the recorded values, and the merged final MIA table (`results/mia_grid_final.csv`) |
| `analysis/mia_followup/` | Rerun of MIA jobs that hit the 300-s grid limit or ended with a process error, under a 3,600-s limit |
| `analysis/metrics_uncertainty/` | Paired, species-pair-clustered summaries of the main grid with cluster-bootstrap intervals; the full-grid accuracy table |
| `analysis/read_length_all_methods/` | Read-length experiment for every method (4X, 50-300 bp and a variable-length condition) |
| `analysis/real_reads_all_methods/` | Every method on the mtDNA-mapped reads of NA07000, scored against the curated 1000 Genomes phase-3 NA07000 sequence (`truth/`), with primary, strict, and no mask, and against rCRS |
| `analysis/runtime_benchmark/` | Controlled run-time and peak-memory benchmark (one job at a time, single thread) |
| `analysis/tie_rule_reruns/` | Reruns of result sets that predated the current tie rule, compared with the v1.3.0 results |
| `analysis/unrefined_block/` | Control-region diagnostics, including the rotation test |
| `figures/` | Scripts that draw the manuscript's data figures (Figs 2-5) |
| `supporting/build_si_tables.py` | Builds S2-S6 Tables (xlsx with a column-description sheet) from the files above |

Generated outputs (`assets/`, `generated/`, `supporting/*.xlsx`) are not
checked in; the commands below recreate them.

## Reproducing the tables and figures from the stored results

Run from this folder with the benchmark environment (`setup_env.sh` in the
repository root; `openpyxl` is needed for the S Tables):

```bash
python analysis/metrics_uncertainty/summarize_main_grid.py
python analysis/metrics_uncertainty/make_table_full_grid.py
python analysis/read_length_all_methods/make_table_read_length.py
python analysis/real_reads_all_methods/make_table_real_reads.py
python analysis/runtime_benchmark/make_table_runtime.py
python figures/make_data_figures.py
python figures/make_fig2_consensus_example.py
(cd supporting && python build_si_tables.py)
```

These steps read only stored CSVs and finish in minutes. Rerunning the
experiments themselves (`run_*.py`, `rerun_external.py`) needs the external
tools (MIA, Pilon, NOVOPlasty, GetOrganelle, BWA, samtools) and the data from
`fetch_data.py` and `extract_real_wgs_reads.sh`; the full reruns took hours.

## Provenance notes

- The experiments were run from an unmodified `git archive` export of v1.3.0
  (`~/plos-reruns/bench-v1.3.0`). The `manifest_*.json` files record that
  location; the scripts here point to `scripts/` in this repository instead.
- MIA's process errors (aborts and segmentation faults) were intermittent; the
  job-level history of every MIA run is kept in `analysis/mia_followup/results/`
  and `analysis/evaluator_n_fix/results/`.
- `verify_unchanged.py` compares F1 values, not sequences: the original runs
  did not retain their outputs, so a matching F1 is consistent with, but does
  not prove, an identical sequence. Repeated reruns of the same job reproduced
  the stored sequence exactly (`results/determinism_check.txt`).
