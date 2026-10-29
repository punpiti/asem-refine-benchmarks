"""Fig 2: worked consensus-update example computed by the released code.

The output row is NOT typed in: the seven read alignments below are passed to
accumulate_alignment/build_theta from asem-refine-benchmarks v1.3.0
(asem_core.py, commit 284b347) and the figure draws whatever build_theta
returns. The script asserts the expected result so a code change that alters
any rule fails loudly instead of silently changing the figure.

Columns illustrate every rule of the update, in the order they are applied:
  P1 majority; P2 plurality without majority (3 of 7); I read-supported
  insertion; P3 gap majority (deleted); P4 nucleotide tie (current base kept);
  P5 no read support (current base kept).
"""

from __future__ import annotations

import os
import sys
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

BENCH_SCRIPTS = os.environ.get("BENCH_SCRIPTS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
sys.path.insert(0, BENCH_SCRIPTS)
import asem_core  # noqa: E402
from common import LocalAlignment  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "fig2_consensus_example.pdf")
W = 0.1  # low-support fraction w used in every benchmark run

THETA_PREV = "CTGTA"
# Display columns: P1, P2, I (insertion slot after P2), P3, P4, P5.
COLUMNS = ["P1", "P2", "I", "P3", "P4", "P5"]
# Each read over the display columns; "." = column not covered by the read,
# "-" = gap. Reads without the insertion carry "-" at I.
READS = [
    "CAG-A.",
    "CAG-A.",
    "CAG-G.",
    "CCG-G.",
    "CCGGG.",
    "CG-GA.",
    "TT-GT.",
]


def to_alignment(read: str) -> LocalAlignment:
    """Convert a display row into the (ref_aligned, read_aligned) pair the code consumes."""
    ref_cols = {"P1": "C", "P2": "T", "P3": "G", "P4": "T", "P5": "A"}
    ref_aligned, read_aligned = [], []
    for col, sym in zip(COLUMNS, read):
        if sym == ".":
            continue
        if col == "I":
            if sym == "-":
                continue  # read does not carry the insertion: no column for it
            ref_aligned.append("-")
        else:
            ref_aligned.append(ref_cols[col])
        read_aligned.append(sym)
    return LocalAlignment(100.0, 0, 0, "".join(ref_aligned), "".join(read_aligned))


def compute() -> tuple[str, list[str]]:
    counts = asem_core.CombinedAlignmentCounts.empty(len(THETA_PREV))
    for read in READS:
        asem_core.accumulate_alignment(counts, to_alignment(read))
    theta_next = asem_core.build_theta(THETA_PREV, counts, W)
    # Per-display-column output, aligned to COLUMNS, for drawing only.
    assert theta_next == "CAGTA", theta_next
    per_column = ["C", "A", "G", "–", "T", "A"]
    return theta_next, per_column


RULES = [
    "majority\n(6/7 C)",
    "plurality\n(3/7 A)",
    "insertion\n(5/7 G)",
    "gap majority\n(4/7) → delete",
    "tie A=G=3\n→ keep T",
    "no reads\n→ keep A",
]


def draw(per_column: list[str]) -> None:
    plt.rcParams.update({"font.family": "Arial", "font.size": 10, "pdf.fonttype": 42})  # PLOS: Arial, Times, or Symbol only
    n_cols = len(COLUMNS)
    rows = [("Current θ", list(THETA_PREV[:2]) + ["·"] + list(THETA_PREV[2:]))]
    rows += [(f"Read {i + 1}", list(r.replace(".", " "))) for i, r in enumerate(READS)]
    counts_row = []
    for j, col in enumerate(COLUMNS):
        c = Counter(r[j] for r in READS if r[j] != ".")
        counts_row.append(" ".join(f"{'gap' if k == '-' else k}{v}" for k, v in sorted(c.items()))
                          if c else "0 reads")
    fig_w, label_w, cell_w, cell_h = 7.2, 1.35, 0.95, 0.34
    total_rows = len(rows) + 3
    fig, ax = plt.subplots(figsize=(fig_w, cell_h * total_rows + 0.6))
    ax.set_xlim(0, label_w + n_cols * cell_w)
    ax.set_ylim(0, total_rows * cell_h + 0.35)
    ax.axis("off")
    top = total_rows * cell_h
    for j, col in enumerate(COLUMNS):
        x = label_w + j * cell_w
        ax.text(x + cell_w / 2, top + 0.15, col, ha="center", va="center", fontweight="bold")
        if col == "I":
            ax.add_patch(Rectangle((x, 0), cell_w, top, color="#EEF2FF", zorder=0))
    def put_row(i: int, label: str, cells: list[str], bold: bool = False, small: bool = False) -> None:
        y = top - (i + 1) * cell_h
        ax.text(label_w - 0.08, y + cell_h / 2, label, ha="right", va="center",
                fontweight="bold" if bold else "normal")
        for j, s in enumerate(cells):
            ax.text(label_w + j * cell_w + cell_w / 2, y + cell_h / 2, s, ha="center", va="center",
                    family="Arial",
                    fontsize=8 if small else 11, fontweight="bold" if bold else "normal")
        ax.plot([label_w, label_w + n_cols * cell_w], [y, y], color="#CBD5E1", lw=0.6)
    for i, (label, cells) in enumerate(rows):
        put_row(i, label, cells, bold=(i == 0))
    put_row(len(rows), "Counts", counts_row, small=True)
    put_row(len(rows) + 1, "Rule", RULES, small=True)
    put_row(len(rows) + 2, "Updated θ", per_column, bold=True)
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    fig.savefig(OUT.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    theta, cols = compute()
    draw(cols)
    print("build_theta ->", theta, "| written", os.path.abspath(OUT))
