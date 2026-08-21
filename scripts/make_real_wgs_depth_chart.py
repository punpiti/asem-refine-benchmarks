"""Line-chart companion for Table 6 (tbl:real-wgs): F1 vs. sequencing depth
on real human WGS reads refined against the chimpanzee mtDNA reference,
matching the table+figure pattern used for tbl:external-tools/Fig. 9 and
tbl:read-length/Fig. 10.

Values are copied verbatim from bmc/thai-sections/03e-real-wgs-table.tex (and
its English counterpart in bmc/sections/03-results.md) -- update DEPTH/F1
below if the table changes, do not let this drift from it.

Usage: python3 make_real_wgs_depth_chart.py
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid")

SCRIPT_DIR = os.path.dirname(__file__)
OUT_PATH = os.path.join(SCRIPT_DIR, "..", "..", "bmc", "assets", "fig_real_wgs_depth.png")

# Matches bmc/thai-sections/03e-real-wgs-table.tex exactly. 0X is the
# initial, unrefined state (chimp reference vs. real human reads); 1-8X are
# ASEM's final-iteration F1.
DEPTHS = [0, 1, 2, 3, 4, 5, 6, 7, 8]
F1 = [0.9370, 0.9544, 0.9635, 0.9698, 0.9740, 0.9762, 0.9777, 0.9783, 0.9791]

COLOR = "#2a78d6"  # dataviz skill: categorical slot 1 -- single series, no legend needed.


def main() -> None:
    fig, ax = plt.subplots(figsize=(7, 5))

    # 0X (initial, unrefined) drawn as an open marker to distinguish it from
    # the six refined ASEM measurements, same convention as the 1-8X point
    # in fig_external_tools_depth.png.
    ax.plot(DEPTHS[1:], F1[1:], marker="o", markersize=8, linewidth=2, color=COLOR)
    ax.plot(DEPTHS[:2], F1[:2], linewidth=2, color=COLOR)
    ax.scatter([0], [F1[0]], s=140, facecolors="none", edgecolors=COLOR, linewidths=1.5, zorder=5)
    ax.annotate("0X: initial\n(unrefined)", xy=(0, F1[0]), xytext=(1.15, F1[0] + 0.006),
                fontsize=8, color="0.35",
                arrowprops=dict(arrowstyle="-", color="0.35", lw=0.8))

    ax.set_xticks(DEPTHS)
    ax.set_xticklabels(["0X", "1X", "2X", "3X", "4X", "5X", "6X", "7X", "8X"])
    ax.set_xlabel("sequencing depth")
    ax.set_ylabel("F1")
    ax.set_ylim(0.930, 0.985)

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150)
    plt.close(fig)
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
