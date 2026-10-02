# Benchmark software: installation and execution

This document describes the ASEM/ASEM-Hybrid benchmark suite: where each
program comes from, how it is installed, how the wrappers call it, which
parameters are used, and where results are written. The information is taken
from the source code, environment history, and package metadata of the
environments used to produce the reported results.

Throughout this document, commands are run from the repository root
(`asem-refine-benchmarks/`) unless stated otherwise, and
`$MAMBA_ROOT_PREFIX/envs/<env>` denotes the location of a micromamba
environment.

## 1. Overview

| Method | Version used | Environment / location |
|---|---:|---|
| ASEM non-recursive | project source | `scripts/baseline_ieee_access.py` |
| ASEM recursive | project source | `scripts/baseline_ojemb.py` |
| SRSC | project source | `scripts/baseline_ecticon.py` |
| ASEM-Hybrid | project source | `scripts/asem_hybrid.py` and the boundary variant |
| NOVOPlasty | 4.3.5 | `$MAMBA_ROOT_PREFIX/envs/genome/` |
| GetOrganelle | 1.7.7.0 | `$MAMBA_ROOT_PREFIX/envs/getorganelle/` |
| MIA | 1.0 | `$MAMBA_ROOT_PREFIX/envs/mia-assembler/` |
| Pilon | 1.24 | `$MAMBA_ROOT_PREFIX/envs/pilon-benchmark/` |
| MITObim | 1.9.1 | `$MAMBA_ROOT_PREFIX/envs/mitobim/`; installable but not included in the results |

The project's own methods are Python source code and need no separate
installation. NOVOPlasty, GetOrganelle, MIA, Pilon, and MITObim are installed
as packages from [Bioconda](https://bioconda.github.io/), with dependencies from
Bioconda and conda-forge, using `micromamba`.

## 2. Data: sources, download, and storage

Data are not software and are therefore not installed into a micromamba
environment. They are downloaded into `data/`, which is kept separate from the
source code. Do not assume the data are present after a fresh clone; run the
data-preparation scripts first.

### 2.1 Reference genomes from NCBI

`scripts/fetch_data.py` downloads nucleotide FASTA records from NCBI Nucleotide
through the NCBI E-utilities, using the accessions defined in the script:

```bash
cd asem-refine-benchmarks/scripts
micromamba run -n genome python fetch_data.py
```

Endpoint:

```text
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
```

Requests use `db=nuccore`, `rettype=fasta`, and `retmode=text`. The script
pauses 0.34 seconds between requests to stay below the unauthenticated NCBI
E-utilities limit of about 3 requests per second. If a destination file already
exists and is non-empty, it is skipped, so the script can be rerun without
downloading everything again.

Default storage locations:

```text
data/phase1_reproduction/       6 mtDNA species for the main simulation grid
data/ieee2021_mtdna_panel/      30 mtDNA species for the divergence panel
data/refs/                      references for the real-WGS validation
```

The main six-species set:

| Species | NCBI accession | Role |
|---|---|---|
| *Homo sapiens* | NC_012920 | Phase 1 and human rCRS evaluation target |
| *Gorilla gorilla* | NC_001645 | Phase 1 same-family comparison |
| *Saimiri boliviensis* | NC_021966 | Phase 1 same-genus comparison |
| *Saimiri sciureus* | NC_012775 | Phase 1 same-genus comparison |
| *Aotus azarai* | NC_021939 | Phase 1 primate reference/target |
| *Varecia variegata* | NC_012773 | Phase 1 and divergent real-read starting reference |

The 30-species set uses the accessions listed in
`docs/ieee2021_mtdna_panel_manifest.md` and in `MANIFESTS` in
`scripts/fetch_data.py`. The genomes are approximately 16.6–17.3 kb each.

Additional references for the real-WGS validation:

| File | Accession / source | Used as |
|---|---|---|
| `data/refs/chrM_rCRS.fasta` | Human NC_012920 | evaluation truth/target |
| `data/refs/chimp_NC_001643.fasta` | Chimpanzee NC_001643 | same-family starting reference |
| `data/refs/chrM_for_cram.fasta` | Copy of NC_012920 with the header renamed to `chrM` | reference for reading the remote CRAM |

### 2.2 Read simulation for the main benchmark

The simulation does not download prebuilt read files. Reads are generated from
the target mtDNA in each job by `simulate_shotgun_reads()`:

- single-end reads, length 150 bp
- coverage 1–8X for the main ASEM/SRSC/Hybrid grid
- NOVOPlasty/GetOrganelle additionally use 10, 15, 20, and 30X
- 3 replicates in the full grid
- seeded with `stable_seed(target_name, depth, replicate)`
- methods with the same species/depth/replicate therefore receive identical
  reads, read for read

The main grid uses every ordered pair of the six species in which the reference
differs from the target: 30 pairs × 8 depths × 3 replicates = 720 jobs per
method.

For the internal methods, reads are generated in memory. External-tool wrappers
write them as temporary FASTA or FASTQ files under `/tmp`, one directory per
job, deleted when the job finishes. No large simulated-read dataset is stored in
`data/`.

### 2.3 Whole-genome sequencing source for the real-WGS validation

The real-read validation starts from human whole-genome sequencing data, but the
benchmark **does not assemble or refine the nuclear whole genome**. The source
is sample NA07000 from the 1000 Genomes Project 30× resequencing collection:

```text
ENA project:  PRJEB31736
ENA run:      ERR3239279
Sample:       NA07000
Source file:  NA07000.final.cram
Remote URL:   https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239279/NA07000.final.cram
Full size:    approximately 18 GB
```

To avoid downloading the whole-genome CRAM, the script uses a samtools remote
indexed range query to fetch only records mapped to `chrM`:

```bash
cd asem-refine-benchmarks

# Build a reference whose header matches @SQ SN:chrM in the CRAM
awk 'NR==1 {print ">chrM"; next} {print}' \
  data/refs/chrM_rCRS.fasta > data/refs/chrM_for_cram.fasta
micromamba run -n genome samtools faidx data/refs/chrM_for_cram.fasta

# Download the CRAI and extract only the chrM slice as BAM
cd scripts
SAMTOOLS_BIN=$MAMBA_ROOT_PREFIX/envs/genome/bin/samtools \
./extract_real_wgs_reads.sh
```

The core commands inside the script are:

```bash
curl -sf "${CRAM_URL}.crai" -o data/NA07000.final.cram.crai

samtools view -b \
  --reference data/refs/chrM_for_cram.fasta \
  -X "$CRAM_URL" data/NA07000.final.cram.crai \
  chrM -o data/NA07000_chrM.bam
```

Resulting files and their locations in the repository:

```text
data/NA07000.final.cram.crai       1,526,808 bytes
data/NA07000_chrM.bam             91,509,533 bytes
data/refs/chrM_for_cram.fasta
data/refs/chrM_for_cram.fasta.fai
```

The reported extraction yielded 2,254,405 records in the chrM BAM before
filtering. The `run_real_wgs*.py` runners then call:

```bash
samtools view -F 3332 data/NA07000_chrM.bam
```

Flag mask 3332 removes unmapped (`0x4`), secondary (`0x100`), duplicate
(`0x400`), and supplementary (`0x800`) records. Only sequences of exactly
150 bp are kept, giving a usable pool of 1,860,886 reads.

Strand note: the SEQ field of mapped reverse-strand records in a BAM is already
in reference orientation, so the runner does not reverse-complement it again.

### 2.4 Storing data outside the repository

The default is `data/` inside the repository. To keep the BAM on another disk,
pass an output directory to the extraction script and point the runner to it
with `GENOME_WGS_DATA_DIR`:

```bash
cd asem-refine-benchmarks/scripts

CHRM_REF_FASTA=/path/to/data/refs/chrM_for_cram.fasta \
SAMTOOLS_BIN=/path/to/samtools \
./extract_real_wgs_reads.sh /path/to/data

GENOME_WGS_DATA_DIR=/path/to/data \
SAMTOOLS_BIN=/path/to/samtools \
micromamba run -n genome python run_real_wgs_benchmark.py
```

The external directory must keep this layout:

```text
/path/to/data/NA07000_chrM.bam
/path/to/data/refs/chrM_rCRS.fasta
/path/to/data/refs/chimp_NC_001643.fasta
```

The hybrid runner still reads the Varecia file from
`data/phase1_reproduction/Varecia_variegata_NC_012773.fasta` inside the
repository.

## 3. Main environment

The repository provides `setup_env.sh` to create the main environment,
`asem-bench`:

```bash
cd asem-refine-benchmarks
bash setup_env.sh
```

Its core commands are:

```bash
micromamba create -n asem-bench -y \
  -c bioconda -c conda-forge \
  --channel-priority flexible \
  "samtools>=1.10" scikit-bio numpy pandas matplotlib seaborn scipy

micromamba run -n asem-bench pip install parasail
```

The reported results were produced with Python from an environment named
`genome`, containing Python 3.14.6, samtools 1.21, scikit-bio 0.7.3,
parasail 1.3.4, and NOVOPlasty 4.3.5. The commands below therefore use
`micromamba run -n genome`. If you create the environment with `setup_env.sh`,
replace this with `-n asem-bench` and install NOVOPlasty into it as well.

## 4. Internal methods

### ASEM non-recursive

```bash
cd asem-refine-benchmarks
ASEM_GRID_WORKERS=14 \
micromamba run -n genome python scripts/run_grid_baseline1.py
```

### ASEM recursive

```bash
cd asem-refine-benchmarks
ASEM_GRID_WORKERS=14 \
micromamba run -n genome python scripts/run_grid_baseline2.py
```

### SRSC

```bash
cd asem-refine-benchmarks
micromamba run -n genome python scripts/run_grid_baseline3.py
```

These three methods use coverage 1–8X, 3 replicates, and 30 ordered pairs, for
720 jobs per method. `ASEM_GRID_WORKERS` sets the number of concurrent jobs
(default 14 for the ASEM runners). ASEM recursive performs at most 6 refinement
iterations.

### ASEM-Hybrid

The runner supports both `legacy` and `boundary` variants. The source default is
still `legacy`, while the reported results use the boundary-recruited variant,
so it must be selected explicitly:

```bash
cd asem-refine-benchmarks
ASEM_HYBRID_VARIANT=boundary \
ASEM_HYBRID_WORKERS=14 \
micromamba run -n genome python scripts/run_grid_baseline_hybrid.py
```

Without `ASEM_HYBRID_VARIANT=boundary`, the results will not match the reported
results.

### Rerun after the tie-handling change

The consensus rule was changed from an alphabetic `argmax()` to keeping the
existing base when nucleotide counts are tied, and tied candidate insertions are
no longer accepted. The full rerun under this rule is driven by a single shell
script:

```bash
cd asem-refine-benchmarks
bash scripts/run_post_tie_rule_grids.sh
```

The script performs these steps automatically, in order:

1. Runs all unit tests.
2. Runs ASEM non-recursive (720 jobs).
3. Runs ASEM recursive (720 jobs).
4. Runs boundary ASEM-Hybrid (720 jobs).
5. Checks that each method has all 720 jobs and builds a before/after comparison.

The default is 14 workers. To change the number of cores:

```bash
ASEM_GRID_WORKERS=10 bash scripts/run_post_tie_rule_grids.sh
```

The runners write each result as soon as its job finishes and read existing job
keys before starting, so the run **can be resumed**: if the process is stopped,
rerun the same command without deleting the CSV, and completed jobs are skipped.
Do not delete or rename `grid_results.post_tie_rule_partial_20260929.csv` while
resuming.

To run without keeping a terminal open:

```bash
cd asem-refine-benchmarks
mkdir -p scripts/results/post_tie_rule_logs
nohup bash scripts/run_post_tie_rule_grids.sh \
  > scripts/results/post_tie_rule_logs/controller.log 2>&1 &
echo $! > scripts/results/post_tie_rule_logs/controller.pid
```

Check progress with:

```bash
bash scripts/check_post_tie_rule_status.sh
```

Per-method logs:

```text
scripts/results/post_tie_rule_logs/baseline1.log
scripts/results/post_tie_rule_logs/baseline2.log
scripts/results/post_tie_rule_logs/hybrid_boundary.log
```

Post-tie-rule results are written to:

```text
scripts/results/baseline1_ieee_access/grid_results.post_tie_rule_partial_20260929.csv
scripts/results/baseline2_ojemb/grid_results.post_tie_rule_partial_20260929.csv
scripts/results/baseline_hybrid_boundary/grid_results.post_tie_rule_partial_20260929.csv
```

When all jobs are complete, `summarize_post_tie_rule.py` creates:

```text
scripts/results/post_tie_rule_summary.json
```

This file reports job/row counts, the jobs whose final result changed, the
deltas in F1/identity/recall/precision, and the summary values used in the
reported tables.

## 5. NOVOPlasty

### Download and installation

Install the `novoplasty` package from Bioconda:

```bash
micromamba install -n genome -y \
  -c bioconda -c conda-forge \
  novoplasty=4.3.5
```

Package source: <https://anaconda.org/bioconda/novoplasty>

Executable:

```text
$MAMBA_ROOT_PREFIX/envs/genome/bin/NOVOPlasty4.3.5.pl
```

### Invocation and parameters

`scripts/ext_novoplasty.py` creates a scratch directory per job, writes the
reads as FASTQ, and calls:

```bash
perl $MAMBA_ROOT_PREFIX/envs/genome/bin/NOVOPlasty4.3.5.pl \
  -c config.txt
```

Configuration:

```text
Type                    = mito
Genome Range            = 12000-22000
K-mer                   = 33
Max memory              = 4 GB
Seed Input              = seed_ref.fasta
Reference sequence      = seed_ref.fasta
Read Length             = 150
Insert size             = 300
Platform                = illumina
Single/Paired           = SE
Combined reads          = reads.fastq
Insert size auto        = yes
Use Quality Scores      = no
```

The wrapper takes the longest sequence among the Contigs, Circularized, Merged,
or Option FASTA files as the output assembly.

### Running the grid

```bash
cd asem-refine-benchmarks
micromamba run -n genome python scripts/run_grid_ext_novoplasty.py
```

- coverage 1–8, 10, 15, 20, and 30X
- 3 replicates per ordered pair/depth
- per-job timeout read from `NOVOPLASTY_TIMEOUT_S` (default 300 seconds)
- up to 14 concurrent jobs, as set in the runner

NOVOPlasty is a de novo organelle assembler. Its failures at 1–8X should be
interpreted as results under ultra-low coverage and this configuration, not as
a general conclusion that NOVOPlasty does not work.

Scoring after strand and circular-origin normalization is described in
Section 9.

## 6. GetOrganelle

### Download and installation

```bash
micromamba create -n getorganelle -y \
  -c bioconda -c conda-forge \
  getorganelle=1.7.7.0
```

Package source: <https://anaconda.org/bioconda/getorganelle>

The environment contains GetOrganelle 1.7.7.0, SPAdes 3.15.5, and
Bowtie2 2.5.4. The program is located at:

```text
$MAMBA_ROOT_PREFIX/envs/getorganelle/bin/get_organelle_from_reads.py
```

### Invocation and parameters

`scripts/ext_getorganelle.py` prepends this environment to `PATH` so that the
matching Bowtie2 and SPAdes are found, and then calls:

```bash
$MAMBA_ROOT_PREFIX/envs/getorganelle/bin/python \
  $MAMBA_ROOT_PREFIX/envs/getorganelle/bin/get_organelle_from_reads.py \
  -u reads.fastq \
  -s seed_ref.fasta \
  -F animal_mt \
  -o out \
  -R 10 \
  -k 21,45,65,85,105 \
  -t 1
```

- `-u`: single-end reads
- `-s`: starting seed
- `-F animal_mt`: animal mitochondrial genome
- `-R 10`: 10 extension rounds
- `-k 21,45,65,85,105`: k-mer set
- `-t 1`: one thread per job, because the grid runs many jobs concurrently

The wrapper takes the longest sequence from `*path_sequence.fasta`.

### Running the grid

```bash
cd asem-refine-benchmarks
GETORGANELLE_ENV_BIN=$MAMBA_ROOT_PREFIX/envs/getorganelle/bin \
micromamba run -n genome python scripts/run_grid_ext_getorganelle.py
```

Coverage 1–8, 10, 15, 20, and 30X; 3 replicates; timeout 600 seconds per job;
up to 14 concurrent jobs.

Scoring after strand and circular-origin normalization is described in
Section 9.

## 7. MIA (Mapping Iterative Assembler)

### Download and installation

```bash
micromamba create -n mia-assembler -y \
  -c conda-forge -c bioconda \
  --channel-priority flexible \
  mapping-iterative-assembler=1.0
```

Package source: <https://anaconda.org/bioconda/mapping-iterative-assembler>

The package file comes from `conda.anaconda.org/bioconda/linux-64/` as
`mapping-iterative-assembler-1.0-h503566f_7.conda`.

Executables:

```text
$MAMBA_ROOT_PREFIX/envs/mia-assembler/bin/mia
$MAMBA_ROOT_PREFIX/envs/mia-assembler/bin/ma
```

### Invocation and parameters

`scripts/ext_mia.py` writes the reference and reads as FASTA and calls:

```bash
mia \
  -r reference.fasta \
  -f reads.fasta \
  -m assembly \
  -c -D -F
```

- `-r`: starting reference
- `-f`: input reads
- `-m`: output root
- `-c`: circular-reference mode
- `-D`: distant-reference mode; keeps low-scoring reads during iterative assembly
- `-F`: write only the final MALN assembly

The wrapper then selects the highest-numbered `assembly.<iteration>` and
converts the MALN file to FASTA:

```bash
ma -M assembly.<final_iteration> -f 5 -I MIA_consensus
```

Note: `ma -f 5` means output format 5 (FASTA); it is not a minimum mapping
quality of 5.

### Running the full grid

MIA was run on the full grid: 30 ordered pairs, depths 1–8X, 3 replicates
(720 jobs):

```bash
cd asem-refine-benchmarks
MIA_ENV_BIN=$MAMBA_ROOT_PREFIX/envs/mia-assembler/bin \
MIA_SCOPE=full \
MIA_TIMEOUT_S=300 \
MIA_PARALLEL_RUNS=6 \
micromamba run -n genome python3 scripts/run_grid_ext_mia.py
```

Results:

```text
scripts/results/ext_mia/grid_results.csv
```

`MIA_TIMEOUT_S=300` must be set, because the reported results use a 300-second
cutoff, while the runner's default is 180 seconds.

### Representative panel (default scope, kept for reference)

The script's default scope is a 12-job representative panel; its results are
kept for reference in `scripts/results/ext_mia_representative/`:

```bash
cd asem-refine-benchmarks
MIA_ENV_BIN=$MAMBA_ROOT_PREFIX/envs/mia-assembler/bin \
MIA_SCOPE=representative \
MIA_TIMEOUT_S=300 \
MIA_PARALLEL_RUNS=6 \
micromamba run -n genome python scripts/run_grid_ext_mia.py
```

The representative panel:

```text
Saimiri_boliviensis -> Saimiri_sciureus    same genus
Homo_sapiens        -> Gorilla_gorilla      same family
Homo_sapiens        -> Varecia_variegata    same order
coverage: 1, 2, 4, and 8X
replicate: 0
```

## 8. Pilon

### Download and installation

```bash
micromamba create -n pilon-benchmark -y \
  -c conda-forge -c bioconda \
  pilon=1.24 bwa=0.7.19 samtools=1.24
```

Package sources:

- <https://anaconda.org/bioconda/pilon>
- <https://anaconda.org/bioconda/bwa>
- <https://anaconda.org/bioconda/samtools>

Paths used by the wrapper:

```text
$MAMBA_ROOT_PREFIX/envs/pilon-benchmark/bin/bwa
$MAMBA_ROOT_PREFIX/envs/pilon-benchmark/bin/samtools
$MAMBA_ROOT_PREFIX/envs/pilon-benchmark/bin/java
$MAMBA_ROOT_PREFIX/envs/pilon-benchmark/share/pilon-1.24-0/pilon.jar
```

### Invocation and parameters

Pilon takes a BAM mapped to the starting reference, so `scripts/ext_pilon.py`
runs:

```bash
bwa index reference.fasta
bwa mem -t 1 reference.fasta reads.fasta > aligned.sam
samtools sort -o aligned.bam aligned.sam
samtools index aligned.bam

java -Xms512m -Xmx1g \
  -jar pilon.jar \
  --genome reference.fasta \
  --unpaired aligned.bam \
  --output pilon \
  --outdir output_directory \
  --fix all
```

- `bwa mem -t 1`: one mapping thread per job
- `--unpaired`: the simulated reads are single-end
- `--fix all`: enables all correction categories supported by Pilon
- `-Xms512m -Xmx1g`: Java heap of 512 MB–1 GB per job

### Running the full grid

Pilon was run on the full grid: 30 ordered pairs, depths 1–8X, 3 replicates
(720 jobs):

```bash
cd asem-refine-benchmarks
PILON_SCOPE=full \
PILON_TIMEOUT_S=300 \
micromamba run -n genome python3 scripts/run_grid_ext_pilon.py
```

Results:

```text
scripts/results/ext_pilon/grid_results.csv
```

### Representative panel (default scope, kept for reference)

The script's default scope is the same 12-job representative panel as MIA
(three representative pairs, coverage 1, 2, 4, and 8X, replicate 0); its
results are kept for reference in `scripts/results/ext_pilon_representative/`:

```bash
cd asem-refine-benchmarks
PILON_PARALLEL_RUNS=6 \
PILON_TIMEOUT_S=180 \
micromamba run -n genome python scripts/run_grid_ext_pilon.py
```

## 9. Strand- and origin-normalized scoring of NOVOPlasty and GetOrganelle

De novo assemblies from NOVOPlasty and GetOrganelle may be reverse-complemented
or start at a different position on the circular genome. Their outputs are
therefore scored after strand and circular-origin normalization.

Because the original wrappers did not keep the assemblies, every originally
successful job was rerun and rescored:

```bash
cd asem-refine-benchmarks
NOVOPLASTY_TIMEOUT_S=3600 \
python3 scripts/rerun_external_successes_normalized.py novoplasty --all-divergences
python3 scripts/rerun_external_successes_normalized.py getorganelle
python3 scripts/summarize_normalized_external.py
```

Outputs:

```text
scripts/results/ext_novoplasty/normalized_successful_rerun.csv
scripts/results/ext_getorganelle/normalized_successful_rerun.csv
scripts/results/ext_novoplasty/assemblies/          retained assemblies
scripts/results/ext_getorganelle/assemblies/        retained assemblies
scripts/results/external_normalized_summary.json    summary
scripts/results/normalized_rerun_20261001/          acceptance records
```

Success/failure status is still taken from each tool's original
`grid_results.csv`; the rerun only rescores the jobs that originally succeeded.

## 10. MITObim: installable but not included in the benchmark

```bash
micromamba create -n mitobim -y \
  -c bioconda -c conda-forge \
  mitobim=1.9.1
```

Package source: <https://anaconda.org/bioconda/mitobim>

The environment contains MITObim 1.9.1, MIRA 4.0.2, Perl 5.22, and Python 2.7.
Executable:

```text
$MAMBA_ROOT_PREFIX/envs/mitobim/bin/MITObim.pl
```

However, the `mirabait` binary bundled with MIRA uses legacy `vsyscall` and
crashes with a segmentation fault on the WSL2 kernel 6.6 used for the
benchmark. The crash occurs even with a one-read smoke test, and MIRA 4.9.6 was
also tested separately. There is therefore no MITObim grid result. This is a
technical exclusion, not a finding that MITObim fails at low coverage.

## 11. Real-WGS benchmark and supplementary experiments

### 11.1 Original real-WGS benchmark

This experiment uses chimpanzee mtDNA `NC_001643` as the starting reference,
human rCRS `NC_012920` as the evaluation target, and subsamples reads without
replacement from the filtered NA07000 pool:

```bash
cd asem-refine-benchmarks
GENOME_WGS_DATA_DIR=<repo>/data \
SAMTOOLS_BIN=$MAMBA_ROOT_PREFIX/envs/genome/bin/samtools \
micromamba run -n genome python scripts/run_real_wgs_benchmark.py
```

Parameters:

- methods: ASEM non-recursive and ASEM recursive
- starting reference: chimpanzee mtDNA NC_001643
- truth: human rCRS NC_012920
- depths: 1–8X
- 3 replicates
- read length: exactly 150 bp
- maximum iterations: 6
- reads per job: `round(depth × 16,569 / 150)`, i.e. about 110 reads at 1X up
  to 884 reads at 8X

Results:

```text
scripts/results/real_wgs_benchmark/grid_results.csv
```

### 11.2 Boundary-Hybrid real-WGS benchmark

This experiment uses the same read subsets to compare two starting references:

1. Chimpanzee NC_001643 → human (same family): checks that Hybrid does not
   degrade cases where ordinary mapping already works well.
2. *Varecia variegata* NC_012773 → human (same order): tests the benefit of
   boundary recruitment when the reference is more distant.

Command matching the reported results:

```bash
cd asem-refine-benchmarks
GENOME_WGS_DATA_DIR=<repo>/data \
SAMTOOLS_BIN=$MAMBA_ROOT_PREFIX/envs/genome/bin/samtools \
ASEM_HYBRID_VARIANT=boundary \
ASEM_REAL_WGS_WORKERS=2 \
micromamba run -n genome python scripts/run_real_wgs_hybrid_benchmark.py
```

Additional parameters:

- ASEM mapping threshold `tau = 0.5`
- `boundary_flank = 150` for the boundary variant
- maximum iterations: 6
- 1 worker inside each algorithm job
- `ASEM_REAL_WGS_WORKERS=2` controls the number of concurrent jobs
- the Pan reference runs ASEM non-recursive, ASEM recursive, and ASEM-Hybrid
- the Varecia reference runs ASEM non-recursive and ASEM-Hybrid

Results:

```text
scripts/results/real_wgs_hybrid_boundary/final_results.csv
```

Do not omit `ASEM_HYBRID_VARIANT=boundary`: the runner default is still
`legacy`, which writes to a different results directory.

### 11.3 Scope of "whole/full genome"

The NA07000 source data are genuine whole-genome sequencing data, but this
benchmark uses only the mitochondrial-mapped read slice and evaluates a
mitochondrial reference of about 16.6 kb. It does not test reconstruction or
refinement of the roughly 3 Gb nuclear genome. The experiment is therefore a
"real-WGS-derived mitochondrial-read validation" (a real-read mtDNA validation),
not a "full-genome benchmark".

Supporting claims at the whole-nuclear-genome level would require a separate
benchmark design, including chromosome/region selection, memory scaling, repeat
handling, structural variation, and an evaluation truth, none of which are in
this repository.

### 11.4 Supplementary same-genus pair (Macaca)

Setting `ASEM_PAIR_SET=macaca_same_genus` makes every grid runner use the
verified same-genus pair *Macaca silenus* (NC_025221) and *Macaca tonkeana*
(NC_025222) from `data/ieee2021_mtdna_panel/`, in both directions. The runners
write to their usual result directories; the published results for this pair
were then moved to:

```text
scripts/results/same_genus_macaca/<method>/
```

For example:

```bash
cd asem-refine-benchmarks
ASEM_PAIR_SET=macaca_same_genus \
micromamba run -n genome python scripts/run_grid_baseline1.py
```

Run this in a separate copy of the repository or with empty result
directories, because the runners append to and resume from existing
`grid_results.csv` files.

### 11.5 Mixed-orientation test

```bash
cd asem-refine-benchmarks
python3 scripts/run_orientation_experiment.py
```

This requires the `asem-refine` package v0.4.1 to be installed. It reruns
Homo_sapiens -> Varecia_variegata at 1, 2, 4, and 8X × 3 replicates for ASEM and
ASEM-Hybrid, both with the grid reads and with each read reverse-complemented
with probability 0.5. Output:

```text
scripts/results/orientation_experiment/results.csv
```

### 11.6 Other supplementary experiments

```bash
cd asem-refine-benchmarks

# Read-length experiment
micromamba run -n genome python scripts/run_read_length_experiment.py

# Computational-cost profiling
micromamba run -n genome python scripts/profile_computational_cost.py

# Supplementary diagnostics
micromamba run -n genome python scripts/run_parameter_sensitivity.py
micromamba run -n genome python scripts/run_circular_boundary_diagnostic.py
micromamba run -n genome python scripts/run_evaluator_endgap_sanity.py
```

## 12. Results and resuming

Runners write CSV files under:

```text
scripts/results/<benchmark-name>/grid_results.csv
```

Examples:

```text
scripts/results/baseline1_ieee_access/grid_results.csv
scripts/results/baseline2_ojemb/grid_results.csv
scripts/results/baseline_hybrid_boundary/grid_results.csv
scripts/results/ext_mia/grid_results.csv
scripts/results/ext_pilon/grid_results.csv
scripts/results/ext_mia_representative/grid_results.csv
scripts/results/ext_pilon_representative/grid_results.csv
```

Many runners support resuming: they read the existing CSV and skip
`(reference, target, depth, replicate)` combinations already present. A repeated
call may therefore find no pending jobs and is not a fresh rerun. For a full
fresh rerun, work in a copy of the repository or move the existing results
elsewhere first. Do not overwrite the raw results without keeping a copy.

The `make_*.py` scripts read the raw CSVs to build figures and tables, for
example:

```bash
micromamba run -n genome python scripts/make_external_tools_depth_chart_english.py
micromamba run -n genome python scripts/make_hybrid_comparison_figure.py
micromamba run -n genome python scripts/make_computational_cost_table.py
micromamba run -n genome python scripts/make_report.py
```

### Manuscript figures

```bash
cd asem-refine-benchmarks
MANUSCRIPT_FIGURE_DIR=<out dir> \
python3 scripts/manuscript_figures/make_vector_figures.py
python3 scripts/manuscript_figures/make_worked_example.py
```

- `make_vector_figures.py` draws every data figure as vector PDF from the
  result CSVs, writing to `MANUSCRIPT_FIGURE_DIR`.
- `make_worked_example.py` draws the worked-example figure.
- `scripts/manuscript_figures/fig1_overview.svg` is the hand-drawn overview
  schematic.

## 13. Snapshot contents

Snapshot `v1.3.0` adds, relative to `v1.2.0`: the full 720-job MIA and Pilon
grids; strand- and origin-normalized rescoring of every successful NOVOPlasty
and GetOrganelle job, with retained assemblies and acceptance records; the
supplementary *Macaca* same-genus pair for all eight methods; the
mixed-orientation test; the scripts that draw the manuscript figures; the
`grid_run_guard.py` helper required by `run_grid_ext_novoplasty.py`; and the
`NOVOPLASTY_TIMEOUT_S` and `ASEM_PAIR_SET` options. Earlier contents (tie-safe
consensus, post-tie-rule grids, parameter sensitivity, circular-boundary
diagnostic, evaluator sanity check) are unchanged. The Hybrid commands specify
`ASEM_HYBRID_VARIANT=boundary`, `ma -f 5` is documented as the FASTA output
format, and the MIA timeout is 300 seconds.

## 14. Verifying the installation

```bash
micromamba env list

$MAMBA_ROOT_PREFIX/envs/genome/bin/NOVOPlasty4.3.5.pl --help
$MAMBA_ROOT_PREFIX/envs/getorganelle/bin/get_organelle_from_reads.py --version
$MAMBA_ROOT_PREFIX/envs/mia-assembler/bin/mia
$MAMBA_ROOT_PREFIX/envs/mia-assembler/bin/ma
$MAMBA_ROOT_PREFIX/envs/pilon-benchmark/bin/java \
  -jar $MAMBA_ROOT_PREFIX/envs/pilon-benchmark/share/pilon-1.24-0/pilon.jar \
  --version
```

Some older programs may lack `--version` or return a non-zero exit status when
asked for help, so check both the executable output and the package metadata
under `<environment>/conda-meta/`.
