#!/usr/bin/env python3
"""PLOS ONE data figures (copy of asem-refine-benchmarks
scripts/manuscript_figures/make_vector_figures.py, adapted).

Changes from the released generator:
  * MIA comes from the N-safe rerun of all 720 jobs
    (analysis/evaluator_n_fix/results/mia_grid_final.csv).
  * The reference-guided depth figure shows F1 and identity for every tier
    instead of an MIA output-rate panel.
  * The read-length figure is drawn from the all-method read-length run.
Only these figures are regenerated; the others in assets/ are unchanged.

Run with: python plos-one/figures/make_data_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colorbar
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os

PLOS = Path(__file__).resolve().parents[1]
BENCH = PLOS.parent / "scripts"
RESULTS = BENCH / "results"
ANALYSIS = PLOS / "analysis"
MIA_FINAL = ANALYSIS / "evaluator_n_fix" / "results" / "mia_grid_final.csv"
READ_LENGTH = ANALYSIS / "read_length_all_methods" / "results" / "read_length_summary.csv"
ASSETS = Path(os.environ.get("MANUSCRIPT_FIGURE_DIR", str(PLOS / "assets")))
ASSETS.mkdir(parents=True, exist_ok=True)
KEYS = ["reference", "target", "depth", "replicate"]



def _vector_rc() -> None:
    # Embed TrueType fonts (Type 42) so text stays editable/searchable text.
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    # Colorbars rasterize their gradient once it has >= n_rasterize steps;
    # keep the colour ramp as vector quads too.
    matplotlib.colorbar.Colorbar.n_rasterize = 10**9


# Print style shared by the manuscript data figures (PLOS: Arial, Times, or
# Symbol at 8-12 pt at the printed size; figures are drawn at <= 7.5 in wide).
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8986", "#e8e7e4"


def _plos_style() -> None:
    plt.rcdefaults()
    _vector_rc()
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8.5,
        "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "axes.labelcolor": INK2,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 3, "ytick.major.size": 3,
    })


def _letter(ax, letter: str, x: float = -0.02) -> None:
    """PLOS-style panel letter: bold capital at the top left of the panel."""
    ax.text(x, 1.04, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            color=INK, ha="right", va="bottom")


def _line_style(name: str, color: str, marker: str, ls: str) -> dict:
    """Filled markers with a white ring; open markers for conditional (success-only) means."""
    conditional = name in ("MIA", "NOVOPlasty", "GetOrganelle")
    return dict(color=color, linestyle=ls, linewidth=1.5, marker=marker, markersize=4.6,
                markerfacecolor="white" if conditional else color,
                markeredgecolor=color if conditional else "white",
                markeredgewidth=1.1 if conditional else 0.8)


def final_iteration(df: pd.DataFrame) -> pd.DataFrame:
    return df.loc[df.groupby(KEYS)["iteration"].idxmax()]


def panel_label(ax, letter: str) -> None:
    """Springer Nature style: bold lowercase part letter at the top left."""
    ax.text(-0.13, 1.06, letter, transform=ax.transAxes, fontsize=14,
            fontweight="bold", ha="left", va="bottom")


def print_fit(fig, size: tuple[float, float], min_font: float = 8.0) -> None:
    """Resize a figure drawn by an external plotting function and raise any
    text below `min_font` so it stays legible at the printed width."""
    from matplotlib.text import Text
    fig.set_size_inches(*size)
    for t in fig.findobj(Text):
        if t.get_text() and t.get_fontsize() < min_font:
            t.set_fontsize(min_font)
    fig.tight_layout()


# ---------------------------------------------------------------------------
# Figures reproduced with the released benchmark plotting functions
# ---------------------------------------------------------------------------

def _seamless(mesh) -> None:
    # Adjacent vector quads are antialiased independently, which leaves
    # hairline seams in PDF viewers; stroking each quad in its own fill
    # colour closes them.
    mesh.set_edgecolor("face")
    mesh.set_linewidth(0.3)


def fig_initial_state_heatmap() -> None:
    """30-species initial-state F1 matrix, phylogeny-clustered.

    Same construction as make_report._clustermap (average linkage on the
    symmetrised F1 distance), applied to results/full_panel.
    """
    import seaborn as sns
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform

    from full_panel_species import list_species

    _vector_rc()
    df = pd.read_csv(RESULTS / "full_panel" / "initial_state_30species.csv")
    order = list(list_species().keys())
    if len(order) != 30 or len(df) != 900:
        raise ValueError(f"expected 30 species / 900 pairs, got {len(order)} / {len(df)}")
    mat = df.pivot_table(index="reference", columns="target", values="f1", aggfunc="mean")
    mat = mat.reindex(index=order, columns=order)
    sym = (mat.values + mat.values.T) / 2.0
    dist = sym.max() - sym
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2.0
    link = linkage(squareform(dist, checks=False), method="average")
    # mako (wide lightness range, no yellow) keeps this narrow-range
    # similarity scale visibly distinct from the 0--1 viridis scale of the
    # cross-method heatmaps, so the two are not read as one scale.
    g = sns.clustermap(mat, row_linkage=link, col_linkage=link, cmap="mako_r",
                       xticklabels=order, yticklabels=order, figsize=(6.4, 6.2),
                       cbar_kws={"label": "initial F1 (no refinement)"}, cbar_pos=(0.90, 0.03, 0.02, 0.14))
    _seamless(g.ax_heatmap.collections[0])
    _seamless(g.ax_cbar.collections[0])
    g.ax_heatmap.set_xlabel("target", fontsize=9)
    g.ax_heatmap.set_ylabel("reference", fontsize=9)
    g.ax_heatmap.tick_params(labelsize=7.5)
    g.ax_cbar.tick_params(labelsize=8)
    g.ax_cbar.yaxis.label.set_size(8.5)
    g.savefig(ASSETS / "fig_initial_state_heatmap.pdf")
    plt.close("all")


DEPTH_COLORS = {4: "#86b6ef", 8: "#1c5cab"}  # validated ordinal blue ramp


def fig_divergence_boxplot() -> None:
    """ASEM with recursion at 4X and 8X: box + every job as a point.

    Same data selection as make_algo_figures._divergence_boxplot (final
    iteration of baseline2_ojemb), redrawn for print: the closer divergence
    classes have only six jobs per depth, so points are shown, not just boxes.
    """
    _vector_rc()
    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.rcParams["font.size"] = 11
    grid = final_iteration(pd.read_csv(RESULTS / "baseline2_ojemb" / "grid_results.csv"))
    sub = grid[grid["depth"].isin([4, 8])]
    metrics = [("f1", "F1", "F1"), ("identity_pct", "Identity (%)", "Identity")]
    rng = np.random.default_rng(7)
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.9))
    for ax, (col, ylabel, title) in zip(axes, metrics):
        for g, div in enumerate(DIVERGENCES):
            for k, depth in enumerate((4, 8)):
                vals = sub[(sub["divergence"] == div) & (sub["depth"] == depth)][col].to_numpy()
                x = g + (k - 0.5) * 0.36
                color = DEPTH_COLORS[depth]
                ax.boxplot(vals, positions=[x], widths=0.28, showfliers=False, patch_artist=True,
                           boxprops=dict(facecolor=color, alpha=0.30, edgecolor=color, linewidth=1.4),
                           medianprops=dict(color=color, linewidth=2.2),
                           whiskerprops=dict(color=color, linewidth=1.2),
                           capprops=dict(color=color, linewidth=1.2))
                jitter = rng.uniform(-0.09, 0.09, size=len(vals))
                ax.scatter(x + jitter, vals, s=13, color=color, edgecolor="white",
                           linewidth=0.4, zorder=3)
        counts = sub[sub["depth"] == 4].groupby("divergence").size()
        ax.set_xticks(range(3), [f"{d.split('-')[1]}\nn = {counts[d]}"
                                 for d in DIVERGENCES])
        ax.set_xlim(-0.6, 2.6)
        ax.set_xlabel("Shared taxonomic rank")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=12)
        panel_label(ax, "ab"[metrics.index((col, ylabel, title))])
        ax.grid(axis="y", color="0.9", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
    handles = [plt.Line2D([], [], marker="s", linestyle="", markersize=10, color=DEPTH_COLORS[d],
                          label=f"{d}X") for d in (4, 8)]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.03), title="Sequencing depth", title_fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92), w_pad=2.5)
    fig.savefig(ASSETS / "fig_divergence_boxplot.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    plt.rcdefaults()


def fig_hybrid_vs_plain() -> None:
    import make_hybrid_comparison_figure as mod

    _vector_rc()
    mod.OUT_PATH = str(ASSETS / "fig_hybrid_vs_plain_depth.pdf")
    orig_savefig = mod.plt.Figure.savefig

    def savefig(fig, *a, **k):
        for ax in fig.axes:
            ax.set_ylabel("Mean F1, same-order pairs")
            ax.set_xlabel("Sequencing depth (X)")
        print_fit(fig, (5.0, 3.5), 8.5)
        return orig_savefig(fig, *a, **k)
    mod.plt.Figure.savefig = savefig
    try:
        mod.main()
    finally:
        mod.plt.Figure.savefig = orig_savefig


def fig_read_length() -> None:
    """Mean F1 by read length (range over replicates as vertical lines), one panel per method."""
    _vector_rc()
    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.rcParams["font.size"] = 10
    df = pd.read_csv(READ_LENGTH)
    lengths = ["50", "75", "100", "150", "250", "300", "variable"]
    x = {l: i for i, l in enumerate(lengths)}
    x["variable"] = len(lengths) - 1 + 0.6
    panels = [("ASEM", "ASEM"), ("ASEM_Hybrid", "ASEM-Hybrid"), ("SRSC", "SRSC"), ("Pilon", "Pilon"), ("MIA", "MIA")]
    tiers = [("same-genus", "same genus", "#4b5563", "o", ":"), ("same-family", "same family", "#2a78d6", "s", "--"),
             ("same-order", "same order", "#eb6834", "^", "-")]
    fig, axes = plt.subplots(1, 5, figsize=(7.4, 2.9), sharey=True)
    for p_i, (ax, (m, title)) in enumerate(zip(axes, panels)):
        for tier, label, color, marker, ls in tiers:
            sub = df[(df["method"] == m) & (df["divergence"] == tier)].set_index("length_condition")
            sub.index = sub.index.astype(str)
            fixed = [l for l in lengths[:-1] if pd.notna(sub.loc[l, "f1_mean"])]
            ax.plot([x[l] for l in fixed], [sub.loc[l, "f1_mean"] for l in fixed], color=color, marker=marker,
                    linestyle=ls, linewidth=1.4, markersize=4, label=label)
            for l in lengths:
                if pd.notna(sub.loc[l, "f1_mean"]):
                    ax.vlines(x[l], sub.loc[l, "f1_min"], sub.loc[l, "f1_max"], color=color, linewidth=1.0)
                    if l == "variable":
                        ax.plot(x[l], sub.loc[l, "f1_mean"], color=color, marker=marker, markersize=4,
                                linestyle="none")
        ax.axvspan(x["variable"] - 0.4, x["variable"] + 0.4, color="0.94", zorder=0)
        ax.set_xticks([x[l] for l in lengths], lengths[:-1] + ["var."], rotation=60, fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.set_ylim(-0.03, 1.03)
        ax.grid(axis="y", color="0.9", linewidth=0.6)
        ax.spines[["top", "right"]].set_visible(False)
        panel_label(ax, "abcde"[p_i])
    axes[0].set_ylabel("Mean F1 (4X)")
    fig.supxlabel("Read length (bp)", fontsize=10, y=0.13)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.12, 1, 1), w_pad=0.6)
    fig.savefig(ASSETS / "fig_read_length_effect.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Reference-guided cross-method figures (include Pilon and MIA)
# ---------------------------------------------------------------------------

# Fixed categorical order (validated palette: CVD-safe adjacent pairs; markers
# and line styles give a second, non-colour encoding).
METHODS = [
    ("ASEM", RESULTS / "baseline1_ieee_access" / "grid_results.csv", True, "#2a78d6", "o", "-"),
    ("ASEM-Hybrid", RESULTS / "baseline_hybrid_boundary" / "grid_results.csv", True, "#eb6834", "s", "-"),
    ("SRSC", RESULTS / "baseline3_ecticon" / "grid_results.csv", False, "#1baf7a", "^", "-"),
    ("Pilon", RESULTS / "ext_pilon" / "grid_results.csv", False, "#eda100", "D", "-"),
    ("MIA", MIA_FINAL, False, "#e87ba4", "v", "--"),
]
DIVERGENCES = ["same-genus", "same-family", "same-order"]


def load_method(path: Path, iterative: bool) -> pd.DataFrame:
    if path == MIA_FINAL:
        df = pd.read_csv(path)
        done = df["final_status"].eq("completed")
        df.loc[~done, ["f1", "identity_pct"]] = np.nan
        if len(df) != 720 or df.duplicated(KEYS).any():
            raise ValueError(f"{path}: expected 720 unique jobs, found {len(df)}")
        return df
    df = pd.read_csv(path)
    if iterative:
        df = final_iteration(df)
    if len(df) != 720 or df.duplicated(KEYS).any():
        raise ValueError(f"{path}: expected 720 unique jobs, found {len(df)}")
    return df


# De novo assemblers: success/failure comes from the original grid; the F1
# of every successful job is replaced by the strand/circular-origin
# normalized rerun. Plotting refuses to proceed if any
# successful job lacks a normalized score, so normalized and unnormalized
# values are never mixed.

# Accepted normalized reruns live in this repository; override only to
# inspect another rerun copy.
DENOVO_RERUN_RESULTS = Path(os.environ.get("DENOVO_RERUN_RESULTS", str(RESULTS)))
DENOVO = [
    ("NOVOPlasty", "ext_novoplasty", "#008300", "P", ":"),
    ("GetOrganelle", "ext_getorganelle", "#6250d6", "X", ":"),
]


def load_denovo(tool_dir: str, allow_preliminary: bool = False) -> tuple[pd.DataFrame, bool]:
    """Return (jobs, preliminary). preliminary=True means the normalized
    rerun is incomplete and *every* success keeps its original score."""
    grid = pd.read_csv(RESULTS / tool_dir / "grid_results.csv")
    if len(grid) != 1080 or grid.duplicated(KEYS).any():
        raise ValueError(f"{tool_dir}: expected 1080 unique jobs, found {len(grid)}")
    ok = grid["error"].fillna("").eq("") & grid["f1"].notna()
    rerun = pd.read_csv(DENOVO_RERUN_RESULTS / tool_dir / "normalized_successful_rerun.csv")
    rerun = rerun[rerun["error"].fillna("").eq("")]
    if rerun.duplicated(KEYS).any():
        raise ValueError(f"{tool_dir}: duplicate keys in normalized rerun")
    normalized = rerun.set_index(KEYS)["f1"]
    out = grid.set_index(KEYS)
    out["f1"] = np.nan
    success_keys = out.index[ok.to_numpy()]
    missing = success_keys.difference(normalized.index)
    if len(missing):
        if not allow_preliminary:
            raise ValueError(f"{tool_dir}: {len(missing)} successful jobs lack a normalized score")
        print(f"[PRELIMINARY] {tool_dir}: {len(missing)}/{len(success_keys)} successes "
              "not yet normalized; using original scores for all successes")
        out.loc[success_keys, "f1"] = grid.set_index(KEYS).loc[success_keys, "f1"].to_numpy()
        return out.reset_index(), True
    out.loc[success_keys, "f1"] = normalized.loc[success_keys].to_numpy()
    return out.reset_index(), False


def stamp(ax_or_fig, text: str = "PRELIMINARY", **kw) -> None:
    kw.setdefault("transform", getattr(ax_or_fig, "transAxes", None) or ax_or_fig.transFigure)
    ax_or_fig.text(0.5, 0.5, text, ha="center", va="center", rotation=30, fontsize=kw.pop("fontsize", 15),
                   color="#c0392b", alpha=0.85, fontweight="bold", zorder=20,
                   bbox=dict(facecolor="white", alpha=0.75, edgecolor="#c0392b"), **kw)


HEATMAP_SPECIES = ["Homo_sapiens", "Gorilla_gorilla", "Saimiri_boliviensis",
                   "Saimiri_sciureus", "Aotus_azarai", "Varecia_variegata"]
HEATMAP_LABELS = ["H. sapiens", "G. gorilla", "S. boliviensis", "S. sciureus",
                  "A. azarai", "V. variegata"]


def _heatmap_matrix(df: pd.DataFrame) -> tuple[np.ndarray, int, int]:
    sub = df[df["depth"] == 8]
    if len(sub) != 90:
        raise ValueError(f"expected 90 jobs at 8X, found {len(sub)}")
    means = sub.groupby(["reference", "target"])["f1"].mean()
    n = len(HEATMAP_SPECIES)
    a = np.full((n, n), np.nan)
    for i, ref in enumerate(HEATMAP_SPECIES):
        for j, tgt in enumerate(HEATMAP_SPECIES):
            if ref == tgt:
                a[i, j] = 1.0
            elif (ref, tgt) in means.index:
                a[i, j] = means.loc[(ref, tgt)]
    return a, int(sub["f1"].notna().sum()), len(sub)


HEATMAP_ABBR = ["Hsa", "Ggo", "Sbo", "Ssc", "Aaz", "Vva"]


def fig_cross_method_heatmaps(panels: list[tuple[str, pd.DataFrame]],
                              preliminary: frozenset = frozenset()) -> None:
    """Matched 8X heatmaps, 2 x 4 grid sized for the printed text width.

    Species are abbreviated (key in the caption) so the text can stay large;
    self cells are left white, no-output cells are gray with an x.
    """
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Rectangle
    _plos_style()
    # Sequential magnitude: one hue (blue), light -> dark. No-output cells are a
    # neutral gray with an x, so they cannot be read as low F1.
    cmap = LinearSegmentedColormap.from_list(
        "blue_seq", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
    cmap.set_bad("#d4d3cf")
    fig, grid = plt.subplots(2, 4, figsize=(7.5, 4.4), facecolor="white")
    axes = list(grid.ravel())
    n = len(HEATMAP_SPECIES)
    centers = np.arange(n) + 0.5
    mesh = None
    for idx, (ax, (name, df)) in enumerate(zip(axes, panels)):
        a, ok, total = _heatmap_matrix(df)
        np.fill_diagonal(a, np.nan)
        mesh = ax.pcolormesh(np.ma.masked_invalid(a), cmap=cmap, vmin=0, vmax=1,
                             edgecolors="white", linewidth=0.8)
        for i in range(n):
            ax.add_patch(Rectangle((i, i), 1, 1, facecolor="white", edgecolor="white", zorder=3))
            for j in range(n):
                if i != j and np.isnan(a[i, j]):
                    ax.text(j + 0.5, i + 0.5, "×", ha="center", va="center",
                            fontsize=9, color=INK2, zorder=4)
        ax.set_xlim(0, n)
        ax.set_ylim(n, 0)
        ax.set_aspect("equal")
        ax.set_title(name, fontsize=9, color=INK, pad=13)
        ax.text(0.5, 1.02, f"{ok}/{total} outputs", transform=ax.transAxes, ha="center", va="bottom",
                fontsize=8, color=INK2)
        ax.set_xticks(centers, HEATMAP_ABBR, rotation=90, fontsize=8)
        ax.set_yticks(centers, HEATMAP_ABBR if idx % 4 == 0 else [], fontsize=8)
        ax.tick_params(length=0, pad=2)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if name in preliminary:
            stamp(ax, fontsize=9)
    # Eighth slot: colour bar and axis key.
    key = axes[7]
    key.axis("off")
    cax = key.inset_axes([0.02, 0.10, 0.10, 0.82])
    cbar = fig.colorbar(mesh, cax=cax)
    cbar.ax.tick_params(labelsize=8, length=2, color=MUTED, labelcolor=INK2)
    cbar.outline.set_visible(False)
    key.text(0.50, 0.92, "Mean F1, 8X", fontsize=8.5, color=INK, transform=key.transAxes, va="top")
    key.text(0.50, 0.70, "Rows: starting\nreference\nColumns: target", fontsize=8, color=INK2,
             transform=key.transAxes, va="top", linespacing=1.3)
    key.text(0.50, 0.24, "× No output", fontsize=8, color=INK2, transform=key.transAxes, va="top")
    fig.tight_layout(w_pad=0.6)
    fig.subplots_adjust(hspace=0.62)
    fig.savefig(ASSETS / "fig_cross_method_heatmaps.pdf", bbox_inches="tight",
                facecolor="white")
    plt.close(fig)


def fig_reference_guided_depth(data: dict[str, pd.DataFrame]) -> None:
    """Mean F1 (top row) and identity (bottom row) by depth for each tier.

    The y-axes start above zero (line chart, not bars) so that the methods
    separate; the ranges are shared within each row. Line ends are labelled
    directly in the same-order column, where the lines separate.
    """
    _plos_style()
    depths = list(range(1, 9))
    npairs = {"same-genus": 2, "same-family": 2, "same-order": 26}
    rows = [("f1", "Mean F1", (0.5, 1.005), [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
            ("identity_pct", "Mean identity (%)", (30, 100.5), [30, 40, 50, 60, 70, 80, 90, 100])]
    fig, axes = plt.subplots(2, 3, figsize=(7.5, 5.2), sharex=True, sharey="row")
    handles = {}
    for r, (metric, ylabel, ylim, yticks) in enumerate(rows):
        for c, div in enumerate(DIVERGENCES):
            ax = axes[r, c]
            ends = []
            for name, _p, _it, color, marker, ls in METHODS:
                df = data[name]
                y = df[df["divergence"] == div].groupby("depth")[metric].mean().reindex(depths)
                if y.min() < ylim[0]:
                    raise ValueError(f"{name} {div} {metric} below the axis floor {ylim[0]}")
                (h,) = ax.plot(depths, y.values, clip_on=False,
                               zorder={"ASEM": 5, "ASEM-Hybrid": 4}.get(name, 3),
                               label="MIA (completed jobs)" if name == "MIA" else name,
                               **_line_style(name, color, marker, ls))
                handles[name] = h
                ends.append((y.values[-1], name))
            if div == "same-order":
                span = ylim[1] - ylim[0]
                placed = []
                for val, name in sorted(ends):
                    pos = max(val, placed[-1] + 0.055 * span) if placed else val
                    placed.append(pos)
                    leader = abs(pos - val) > 0.01 * span
                    ax.annotate(name, xy=(8, val), xytext=(8.45, pos), fontsize=7.5, color=INK2,
                                va="center", annotation_clip=False,
                                arrowprops=dict(arrowstyle="-", color="#b8b7b3", lw=0.7, shrinkA=1, shrinkB=2)
                                if leader else None)
            ax.set_ylim(*ylim)
            ax.set_yticks(yticks)
            ax.set_xticks(depths)
            ax.set_xlim(0.6, 8.4)
            ax.grid(axis="y", color=GRID, linewidth=0.6)
            ax.set_axisbelow(True)
            ax.spines[["top", "right"]].set_visible(False)
            if r == 0:
                ax.set_title(f"{div.replace('-', ' ').capitalize()} ({npairs[div]} pairs)", color=INK, pad=6)
            else:
                ax.set_xlabel("Mitochondrial depth (X)")
            if c == 0:
                ax.set_ylabel(ylabel)
            _letter(ax, "ABCDEF"[3 * r + c], x=-0.24 if c == 0 else -0.02)
    fig.legend(handles=list(handles.values()), loc="lower center", ncol=5, frameon=False,
               bbox_to_anchor=(0.5, 0.0), handlelength=2.4, columnspacing=1.6)
    fig.tight_layout(rect=(0, 0.05, 0.93, 1), w_pad=1.0, h_pad=1.4)
    fig.savefig(ASSETS / "fig_reference_guided_depth.pdf", facecolor="white")
    plt.close(fig)


DEPTH_AXIS = [1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 30]


def fig_all_methods_depth(data: dict[str, pd.DataFrame],
                          denovo: dict[str, pd.DataFrame],
                          preliminary: frozenset = frozenset()) -> None:
    """Output rate and conditional F1 by depth for all seven methods.

    Reference-guided methods were run at 1--8X only; the de novo assemblers
    additionally at 10--30X. Depths are placed at evenly spaced categorical
    positions with a visible break after 8X.
    """
    _plos_style()
    # Wider spacing after the break keeps the 10/15/20/30 labels apart.
    pos = {d: (i if d <= 8 else 8 + 0.9 + 1.5 * (i - 8)) for i, d in enumerate(DEPTH_AXIS)}
    series = [(n, data[n], c, m, ls) for n, _p, _it, c, m, ls in METHODS]
    series += [(n, denovo[n], c, m, ls) for n, _d, c, m, ls in DENOVO]
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.9))
    for name, df, color, marker, ls in series:
        depths = [d for d in DEPTH_AXIS if d in set(df["depth"])]
        g = df.groupby("depth")["f1"]
        rate = (100 * g.apply(lambda v: v.notna().mean())).reindex(depths)
        cond = g.mean().reindex(depths)
        x = [pos[d] for d in depths]
        style = dict(_line_style(name, color, marker, ls), clip_on=False)
        axes[0].plot(x, rate.values, label=name, **style)
        axes[1].plot(x, cond.values, label=name, **style)
    for ax in axes:
        ax.set_xticks([pos[d] for d in DEPTH_AXIS], [str(d) for d in DEPTH_AXIS])
        ax.axvline((pos[8] + pos[10]) / 2, color=MUTED, linewidth=0.7, linestyle=(0, (2, 2)))
        ax.set_xlabel("Mitochondrial depth (X)")
        ax.grid(axis="y", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylim(-3 if ax is axes[0] else -0.03, 103 if ax is axes[0] else 1.03)
        ax.text((pos[10] + pos[30]) / 2, 1.0, "de novo only",
                transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=8, color=INK2, style="italic")
    axes[0].set_title("Jobs with scoreable output", color=INK, pad=14)
    _letter(axes[0], "A", x=-0.14)
    axes[0].set_ylabel("Jobs with output (%)")
    axes[1].set_title("Mean F1 among jobs with output", color=INK, pad=14)
    _letter(axes[1], "B", x=-0.14)
    axes[1].set_ylabel("Mean F1")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=7, frameon=False,
               bbox_to_anchor=(0.5, 0.0), handlelength=2.2, columnspacing=1.2)
    fig.tight_layout(rect=(0, 0.07, 1, 1), w_pad=2.5)
    if preliminary:
        for ax in axes:
            stamp(ax, "PRELIMINARY\n" + ", ".join(sorted(preliminary)) + " not normalized",
                  fontsize=12)
    fig.savefig(ASSETS / "fig_external_tools_depth.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    fig_read_length()
    plt.rcdefaults()
    data = {name: load_method(path, it) for name, path, it, *_ in METHODS}
    fig_reference_guided_depth(data)
    plt.rcdefaults()
    loaded = {name: load_denovo(d) for name, d, *_ in DENOVO}
    denovo = {n: df for n, (df, _) in loaded.items()}
    fig_cross_method_heatmaps([(n, data[n]) for n, *_ in METHODS] + [(n, denovo[n]) for n, *_ in DENOVO])
    plt.rcdefaults()
    fig_all_methods_depth(data, denovo)
    for name in ("fig_read_length_effect", "fig_reference_guided_depth", "fig_cross_method_heatmaps",
                 "fig_external_tools_depth"):
        print(ASSETS / f"{name}.pdf")


if __name__ == "__main__":
    main()
