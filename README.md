# asem-refine-benchmarks

Current frozen benchmark release: **v1.2.0**.

Reproduction package for the benchmarks reported in the ASEM/ASEM-Hybrid
BMC Genomics paper: the 720-job Phase-1 simulation grid, the
external-tool comparison, the read-length experiment, the computational-cost
profiling, parameter and circular-boundary diagnostics, MIA/Pilon comparisons,
and the real-WGS (NA07000) validation.

This is the exact research code and raw results that produced the numbers
in the paper -- not the `asem-refine` package (https://github.com/punpiti/asem-refine),
which is a separately maintained, independent reimplementation meant for
general use rather than as a frozen benchmark harness. Where the two
differ, this repository is the one that generated the reported tables and
figures.

Version 1.2.0 contains the insertion-aware combined-alignment consensus used
by the revision manuscript and complete result CSVs regenerated with that
implementation. Its consensus update removes a sufficiently supported column
only under a strict gap majority. A unique nucleotide plurality is called;
tied current-reference columns retain their previous base, while tied candidate
insertions are omitted. Read-supported candidate insertion columns can be
retained between reference positions.

## Layout

- `scripts/` -- all benchmark/grid-runner code, the ASEM/ASEM-Hybrid/SRSC
  implementations as actually run, external-tool wrapper scripts, and the
  plotting/table-generation scripts that turned the raw CSVs into the
  paper's figures and tables.
- `scripts/results/` -- the raw result CSVs each grid run actually
  produced. These are the numbers the paper's tables and figures were
  built from; nothing here is synthetic or example data.
- `data/` -- not checked in. Populated by `fetch_data.py` (public NCBI
  reference sequences, downloaded by accession) and, for the real-WGS
  validation only, `extract_real_wgs_reads.sh` (a public ENA read pool).
  See "Getting the data" below.

## Setup

```bash
bash setup_env.sh          # creates a micromamba env with all Python/bio dependencies
micromamba activate asem-bench
```

## Getting the data

Reference sequences are pre-existing public NCBI GenBank/RefSeq records;
this repo does not embed copies of them, only fetches them by accession
into `data/` (gitignored):

```bash
cd scripts
python3 fetch_data.py
```

This populates `data/phase1_reproduction/` (six-species Phase-1 panel),
`data/ieee2021_mtdna_panel/` (30-species divergence-landscape panel), and
`data/refs/` (human rCRS and chimpanzee references for the real-WGS
validation).

The real-WGS validation (Section 3.8) additionally needs the actual read
pool, which is a public human sequencing dataset rather than a small
reference sequence:

```bash
./extract_real_wgs_reads.sh
```

See that script's header comment for what it does and why it only
downloads the mtDNA-mapped slice (a remote-indexed range query), not the
full ~18GB whole-genome CRAM.

## Running the benchmarks

Each `run_grid_*.py` / `run_real_wgs*.py` / `run_read_length_experiment.py`
script is self-contained and writes its output CSV(s) under
`scripts/results/<name>/`, matching what's already checked in there. For
example:

```bash
cd scripts
python3 run_grid_baseline1.py   # non-recursive ASEM, full 720-job grid
python3 run_grid_baseline2.py   # recursive ASEM, full 720-job grid
python3 run_grid_baseline3.py   # SRSC, same grid
ASEM_HYBRID_VARIANT=boundary python3 run_grid_baseline_hybrid.py
```

The complete tie-rule verification can be resumed with one command from the
repository root:

```bash
bash scripts/run_post_tie_rule_grids.sh
bash scripts/check_post_tie_rule_status.sh
```

External-tool comparators (`run_grid_ext_novoplasty.py`,
`run_grid_ext_getorganelle.py`) additionally need NOVOPlasty/GetOrganelle
installed -- see `setup_env.sh`'s optional section.

MIA and Pilon default to the historical 12-job representative panel. Their
shared full-grid runs use the same 30 ordered pairs, depths 1--8X, three
replicates, and deterministic read seeds as the internal methods:

```bash
MIA_SCOPE=full MIA_TIMEOUT_S=300 python3 scripts/run_grid_ext_mia.py
PILON_SCOPE=full PILON_TIMEOUT_S=300 python3 scripts/run_grid_ext_pilon.py
```

Both runners append one row after each completed job and skip keys already in
their result CSV, so interrupted runs can be resumed with the same command.

NOVOPlasty and GetOrganelle outputs are scored after strand and
circular-origin normalization. The original wrappers did not retain
assemblies, so every originally successful job was rerun and rescored
(`normalized_successful_rerun.csv`, assemblies retained under `assemblies/`);
success/failure status still comes from the original `grid_results.csv`.
Some NOVOPlasty successes need more than the first-pass 300-s limit:

```bash
NOVOPLASTY_TIMEOUT_S=3600 python3 scripts/rerun_external_successes_normalized.py novoplasty --all-divergences
python3 scripts/rerun_external_successes_normalized.py getorganelle
python3 scripts/summarize_normalized_external.py
```

Acceptance records for the run used in the paper are in
`scripts/results/normalized_rerun_20261001/`.

The `make_*.py` scripts regenerate the paper's figures/tables from the
CSVs in `scripts/results/`.

## Notes

- Hardcoded machine-specific paths from the original research environment
  (tool binary locations, a personal data directory) have been replaced
  with environment variables that fall back to sensible defaults relative
  to this repo (`SAMTOOLS_BIN`, `NOVOPLASTY_PL`, `GETORGANELLE_ENV_BIN`,
  `GENOME_WGS_DATA_DIR`, `GENOME_ENV_PYTHON`) -- set them if your tools
  live somewhere else or aren't on `PATH`.
- `scripts/asem_hybrid.py` and `scripts/asem_hybrid_boundary.py` are the
  exact algorithm implementations used to generate the reported results,
  kept here even though `asem-refine` now also implements the
  boundary-recruited variant, for byte-for-byte reproducibility of this
  paper's specific numbers.
- License: MIT (see `LICENSE`).
