"""Re-score the stored main-grid NOVOPlasty assemblies that contain ambiguous bases.

GetOrganelle's 198 stored assemblies contain only A/C/G/T, and N-free
sequences are scored identically by both evaluators, so only these files can
change. Prints recorded, original-evaluator and N-safe F1 for each.
"""

from __future__ import annotations

import csv
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import nsafe  # noqa: E402
import common  # noqa: E402
from common import read_fasta  # noqa: E402
from phase1_species import fasta_path  # noqa: E402

DIR = os.path.join(nsafe.BENCH_SCRIPTS, "results", "ext_novoplasty")


def main() -> None:
    recorded = {}
    for path in glob.glob(os.path.join(DIR, "*.csv")):
        for r in csv.DictReader(open(path, newline="")):
            if "f1" in r:
                recorded[(r["reference"], r["target"], int(r["depth"]), int(r["replicate"]))] = r["f1"]
    original = common.pair_align_nucl
    for path in sorted(glob.glob(os.path.join(DIR, "assemblies", "*.fasta"))):
        theta = read_fasta(path)[1]
        if set(theta) <= set("ACGT"):
            continue
        ref, target, d, r = os.path.basename(path)[:-6].split("__")
        key = (ref, target, int(d[1:]), int(r[1:]))
        t = read_fasta(fasta_path(target))[1]
        common.pair_align_nucl = original
        old = common.evaluate_theta_vs_target(theta, t, normalize_strand=True, normalize_circular_origin=True)["f1"]
        nsafe.install()
        new = common.evaluate_theta_vs_target(theta, t, normalize_strand=True, normalize_circular_origin=True)["f1"]
        n_amb = sum(c not in "ACGT" for c in theta)
        print(f"{key} len={len(theta)} ambiguous={n_amb} recorded={recorded.get(key)} old={old:.6f} fixed={new:.6f}")


if __name__ == "__main__":
    main()
