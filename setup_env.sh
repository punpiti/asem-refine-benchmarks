#!/usr/bin/env bash
# One-shot environment setup for this benchmark suite.
# Usage: bash setup_env.sh
# Safe to re-run: micromamba create -y overwrites an existing env of the same name.

set -euo pipefail

command -v micromamba >/dev/null 2>&1 || {
  echo "micromamba not found. Install it first: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html" >&2
  exit 1
}

echo "== Creating/refreshing the 'asem-bench' env (samtools, scikit-bio, numpy, pandas, matplotlib, seaborn, scipy) =="
micromamba create -n asem-bench -y -c bioconda -c conda-forge --channel-priority flexible \
  "samtools>=1.10" scikit-bio numpy pandas matplotlib seaborn scipy
micromamba run -n asem-bench pip install parasail

echo
echo "== Optional: external de novo assemblers, only needed to reproduce Section 2.8 =="
echo "  NOVOPlasty: micromamba create -n asem-bench -c bioconda novoplasty  (adds it to the same env)"
echo "  GetOrganelle: micromamba create -n getorganelle -c bioconda getorganelle  (separate env; set"
echo "                GETORGANELLE_ENV_BIN to its bin/ dir before running run_grid_ext_getorganelle.py)"
echo
echo "'asem-bench' env ready. Known caveats:"
echo "  - EMBOSS/needle is not used here; alignment is parasail Smith-Waterman throughout."
echo "  - NOVOPlasty: invoke as 'perl NOVOPlasty4.3.5.pl -c config.txt', not the .pl directly."
echo "  - MITObim could not be run on a modern kernel when this suite was built (mirabait/glibc/"
echo "    vsyscall incompatibility) -- see the paper's Section 2.8 for detail. Not included here."
