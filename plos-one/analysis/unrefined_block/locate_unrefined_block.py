"""Locate target bases that ASEM leaves unrecovered (u_t) on the same-family pair.

Reruns ASEM (v1.3.0 benchmark code, main-grid seed) for Homo<->Gorilla at 8X,
aligns the output to the target with the published evaluator alignment, and
reports runs of target positions aligned to gaps (missing from the output),
with the rCRS/target feature they fall in.
"""
import os, sys
BENCH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "scripts")
sys.path.insert(0, BENCH)
from common import read_fasta, simulate_shotgun_reads, stable_seed
from baseline_ieee_access import run_asem
from phase1_species import fasta_path
from skbio.alignment import pair_align_nucl
from skbio.sequence import DNA

def runs(ref, tgt, depth=8, rep=0):
    s = read_fasta(fasta_path(ref))[1]; t = read_fasta(fasta_path(tgt))[1]
    reads = simulate_shotgun_reads(t, float(depth), 150, stable_seed(tgt, depth, rep))
    theta = run_asem(s, reads).theta
    a, b = (str(x) for x in pair_align_nucl(DNA(theta), DNA(t), mode="global").paths[0].to_aligned([DNA(theta), DNA(t)]))
    out, tpos, start, extra = [], 0, None, 0
    for x, y in zip(a, b):
        if y != "-":
            tpos += 1
            if x == "-":
                start = tpos if start is None else start
            elif start is not None:
                out.append((start, tpos - 1)); start = None
        else:
            extra += 1
    if start is not None:
        out.append((start, tpos))
    return len(theta), len(t), extra, [(a_, b_) for a_, b_ in out if b_ - a_ + 1 >= 10]

for ref, tgt in [("Gorilla_gorilla", "Homo_sapiens"), ("Homo_sapiens", "Gorilla_gorilla")]:
    lt, ltgt, extra, blocks = runs(ref, tgt)
    print(f"{ref}->{tgt}: theta {lt} target {ltgt} extra-in-theta {extra}; missing runs >=10 bp (target coords):",
          [(a, b, b - a + 1) for a, b in blocks])
