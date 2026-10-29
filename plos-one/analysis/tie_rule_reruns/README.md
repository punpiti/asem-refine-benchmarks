# Tie-rule reruns of v1.3.0 result sets

Six released result sets of asem-refine-benchmarks v1.3.0 predate, or cannot be
shown to postdate, the 2026-09-30 tie rule. They were regenerated on 2026-10-28
from an unmodified `git archive` export of v1.3.0 (`~/plos-reruns/bench-v1.3.0`,
`rerun_tie_rule.sh` and `rerun_real_wgs.sh`). `compare_tie_rule_reruns.py` matches
rows on their key columns and writes `results/comparison.txt`; the rerun CSVs are
copied to `results/`.

Outcome (2026-10-28):

| Set | Matched rows | F1 changed | max abs dF1 | Used in PLOS manuscript |
|---|---|---|---|---|
| read_length_experiment | 189 | 3 | 7.1e-05 | superseded by `read_length_all_methods` (current code) |
| parameter_sensitivity | 396 | 0 | 0 | yes; released values valid |
| circular_boundary | 240 | identical | 0 | yes; released values valid |
| real_wgs_benchmark | 172 | 171 | 0.0011 | no; superseded by `real_reads_all_methods` |
| real_wgs_hybrid (legacy variant) | 120 | 120 | 0.015 | no |
| real_wgs_hybrid_boundary | 120 | 101 | 0.00055 | no; superseded by `real_reads_all_methods` |

The manuscript's empirical-read and read-length results come from the
all-method runs, which use the current (tie-rule) code throughout, so no
manuscript value depends on a pre-tie-rule result set.
