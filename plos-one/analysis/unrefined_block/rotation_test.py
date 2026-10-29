"""Is ASEM's unrefined control region caused by divergence or by the linear origin?

Rotate reference and target by the same offset so the control region sits in
the middle of both linear strings, simulate reads from the rotated target with
the main-grid seed rule, rerun ASEM, and compare unrecovered target bases (u_t)
with the unrotated run.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "scripts"))
from common import evaluate_theta_vs_target, read_fasta, simulate_shotgun_reads, stable_seed
from baseline_ieee_access import run_asem
from phase1_species import fasta_path

OFFSET = 8000  # moves the control region (around the origin) to the middle
for ref, tgt in [("Gorilla_gorilla", "Homo_sapiens"), ("Homo_sapiens", "Gorilla_gorilla")]:
    s0 = read_fasta(fasta_path(ref))[1]; t0 = read_fasta(fasta_path(tgt))[1]
    for label, s, t in [("original", s0, t0), ("rotated", s0[OFFSET:] + s0[:OFFSET], t0[OFFSET:] + t0[:OFFSET])]:
        for rep in range(3):
            reads = simulate_shotgun_reads(t, 8.0, 150, stable_seed(tgt, 8, rep))
            m = evaluate_theta_vs_target(run_asem(s, reads).theta, t)
            print(f"{ref}->{tgt} {label} rep{rep}: F1={m['f1']:.4f} identity={m['identity_pct']:.2f} "
                  f"mismatch={m['w_theta']} missing_target={m['u_t']}", flush=True)
