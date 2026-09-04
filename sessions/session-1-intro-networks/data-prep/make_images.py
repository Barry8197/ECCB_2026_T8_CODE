"""
Generate the explanatory diagrams used in the Session 1 notebooks.

Developer script, like `build_kg_data.py` - participants never run it. Rerun it
if the colour scheme in `s1_helpers.TYPE_COLOURS` changes, so that the diagrams
keep matching the plots the notebooks produce.

    python make_images.py

Outputs to ../images/:
    kg_schema.png            the three node types and how they connect
    curated_vs_computed.png  knowledge graph vs co-expression network
    bipartite_projection.png collapsing genes-and-diseases onto diseases alone
    ontology_climb.png       inheriting an ICD-10 code from an ancestor

`vertex_types.png` in that folder is not generated here - it is reused from
earlier gene co-expression network teaching material.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, str(Path(__file__).parent.parent))
from s1_helpers import TYPE_COLOURS  # noqa: E402

OUT = Path(__file__).parent.parent / "images"
GREY = "#666666"
DPI = 160


def box(ax, xy, w, h, label, colour, fontsize=11, text_colour="white"):
    ax.add_patch(FancyBboxPatch(
        (xy[0] - w / 2, xy[1] - h / 2), w, h,
        boxstyle="round,pad=0.02", facecolor=colour, edgecolor="none"))
    ax.text(xy[0], xy[1], label, ha="center", va="center",
            fontsize=fontsize, color=text_colour, weight="bold")


def arrow(ax, start, end, label="", colour=GREY, style="-|>", offset=0.0,
          fontsize=9, rad=0.0):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle=style, mutation_scale=14, color=colour,
        linewidth=1.6, connectionstyle=f"arc3,rad={rad}"))
    if label:
        mx, my = (start[0] + end[0]) / 2, (start[1] + end[1]) / 2
        ax.text(mx, my + offset, label, ha="center", va="center", fontsize=fontsize,
                color=colour, style="italic",
                bbox=dict(facecolor="white", edgecolor="none", pad=1.5))


def finish(fig, ax, path, xlim, ylim, title=None):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=13, weight="bold", pad=12)
    fig.tight_layout()
    fig.savefig(OUT / path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {path}")


# --------------------------------------------------------------------------

def kg_schema():
    """The schema: what kinds of node exist and which edges join them."""
    fig, ax = plt.subplots(figsize=(9, 4.2))

    box(ax, (1.5, 2.4), 2.0, 0.75, "gene", TYPE_COLOURS["gene"], 13)
    box(ax, (5.0, 2.4), 2.0, 0.75, "disease", TYPE_COLOURS["disease"], 13)
    box(ax, (8.5, 2.4), 2.0, 0.75, "icd10", TYPE_COLOURS["icd10"], 13)

    arrow(ax, (2.55, 2.4), (3.95, 2.4), "associated_with", offset=0.28)
    arrow(ax, (6.05, 2.4), (7.45, 2.4), "maps_to", offset=0.28)
    # is_a is a disease -> disease edge, so it loops back on the same box
    arrow(ax, (4.6, 2.85), (5.4, 2.85), colour=GREY, rad=-1.6)
    ax.text(5.0, 3.72, "is_a  (81 edges)", ha="center", fontsize=9, color=GREY,
            style="italic")

    ax.text(1.5, 1.75, "760 nodes", ha="center", fontsize=9, color=GREY)
    ax.text(5.0, 1.75, "90 nodes", ha="center", fontsize=9, color=GREY)
    ax.text(8.5, 1.75, "31 nodes", ha="center", fontsize=9, color=GREY)

    ax.text(3.25, 2.05, "1,656 edges", ha="center", fontsize=8, color=GREY)
    ax.text(6.75, 2.05, "31 edges", ha="center", fontsize=8, color=GREY)

    # A concrete instance of each, underneath the abstract schema
    ax.plot([0.3, 9.7], [1.25, 1.25], color="#DDDDDD", linewidth=1)
    ax.text(0.3, 0.95, "for example:", fontsize=9, color=GREY, style="italic")
    box(ax, (1.5, 0.42), 2.0, 0.6, "BRCA1", TYPE_COLOURS["gene"], 11)
    box(ax, (5.0, 0.42), 2.3, 0.6, "breast cancer", TYPE_COLOURS["disease"], 11)
    box(ax, (8.5, 0.42), 2.0, 0.6, "C50", TYPE_COLOURS["icd10"], 11)
    arrow(ax, (2.55, 0.42), (3.8, 0.42), "score 0.92", offset=0.24, fontsize=8)
    arrow(ax, (6.2, 0.42), (7.45, 0.42))

    finish(fig, ax, "kg_schema.png", (0, 10), (0, 4.1),
           "The schema: three node types, three edge types")


def curated_vs_computed():
    """Why a knowledge graph is a different object from a co-expression network."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

    # Left: co-expression. One node type, edges computed from data.
    ax = axes[0]
    coords = [(1.2, 3.2), (2.8, 3.6), (2.2, 2.0), (3.6, 2.4), (1.0, 1.6), (3.2, 1.0)]
    for a in range(len(coords)):
        for b in range(a + 1, len(coords)):
            if (a + b) % 2 == 0:
                ax.plot(*zip(coords[a], coords[b]), color="#CCCCCC", linewidth=1.4, zorder=1)
    for x, y in coords:
        ax.scatter(x, y, s=520, color=TYPE_COLOURS["gene"], zorder=2, edgecolors="white")
    ax.text(2.3, 0.25,
            "edges are $\\bf{computed}$\ncorrelation above a threshold,\ncalculated from our data",
            ha="center", fontsize=9.5, color=GREY)
    ax.set_title("Computed network\n(gene co-expression, Part 1 §7)",
                 fontsize=11.5, weight="bold")
    ax.set_xlim(0.2, 4.4); ax.set_ylim(-0.3, 4.3); ax.axis("off")

    # Right: knowledge graph. Several node types, edges asserted by a curator.
    ax = axes[1]
    genes = [(1.0, 3.4), (1.0, 2.5), (1.0, 1.6)]
    diseases = [(2.6, 3.0), (2.6, 1.9)]
    icd = [(4.1, 2.45)]
    for g in genes:
        for d in diseases:
            ax.plot(*zip(g, d), color="#CCCCCC", linewidth=1.4, zorder=1)
    ax.plot(*zip(diseases[0], diseases[1]), color="#AAAAAA",
            linewidth=1.8, linestyle="--", zorder=1)
    ax.plot(*zip(diseases[0], icd[0]), color="#CCCCCC", linewidth=1.4, zorder=1)
    for pts, key in ((genes, "gene"), (diseases, "disease"), (icd, "icd10")):
        for x, y in pts:
            ax.scatter(x, y, s=520, color=TYPE_COLOURS[key], zorder=2, edgecolors="white")
    ax.text(2.6, 0.25,
            "edges are $\\bf{asserted}$\nsomeone read the evidence\nand recorded them",
            ha="center", fontsize=9.5, color=GREY)
    ax.set_title("Knowledge graph\n(Open Targets + MONDO, this session)",
                 fontsize=11.5, weight="bold")
    ax.set_xlim(0.2, 4.9); ax.set_ylim(-0.3, 4.3); ax.axis("off")
    ax.legend(handles=[mpatches.Patch(color=TYPE_COLOURS[k], label=k)
                       for k in ("gene", "disease", "icd10")],
              loc="upper right", fontsize=8.5, frameon=False)

    fig.tight_layout()
    fig.savefig(OUT / "curated_vs_computed.png", dpi=DPI, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print("  wrote curated_vs_computed.png")


def bipartite_projection():
    """Collapsing a two-mode graph onto one of its modes."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    gene_y = [3.4, 2.75, 2.1, 1.45, 0.8]
    labels = ["BRCA1", "BRCA2", "TP53", "PTEN", "EGFR"]
    dis = {"breast": (3.2, 2.9), "ovarian": (3.2, 1.3)}
    links = {"breast": [0, 1, 2, 4], "ovarian": [0, 1, 2, 3]}

    ax = axes[0]
    for name, (dx, dy) in dis.items():
        for i in links[name]:
            ax.plot([1.3, dx], [gene_y[i], dy], color="#CCCCCC", linewidth=1.3, zorder=1)
    for y, label in zip(gene_y, labels):
        ax.scatter(1.3, y, s=430, color=TYPE_COLOURS["gene"], zorder=2, edgecolors="white")
        ax.text(0.95, y, label, ha="right", va="center", fontsize=9)
    for name, (dx, dy) in dis.items():
        ax.scatter(dx, dy, s=560, color=TYPE_COLOURS["disease"], zorder=2, edgecolors="white")
        ax.text(dx + 0.28, dy, name, ha="left", va="center", fontsize=10, weight="bold")
    ax.set_title("Bipartite: genes ── diseases", fontsize=11.5, weight="bold")
    ax.set_xlim(-0.2, 4.6); ax.set_ylim(0.2, 4.0); ax.axis("off")

    ax = axes[1]
    ax.plot([1.4, 3.1], [2.6, 2.6], color=TYPE_COLOURS["disease"], linewidth=4, zorder=1)
    for x, name in ((1.4, "breast"), (3.1, "ovarian")):
        ax.scatter(x, 2.6, s=560, color=TYPE_COLOURS["disease"], zorder=2, edgecolors="white")
        ax.text(x, 2.15, name, ha="center", fontsize=10, weight="bold")
    ax.text(2.25, 2.95, "3 shared genes\nJaccard 0.60", ha="center", fontsize=9.5,
            color=GREY)
    ax.text(2.25, 1.15,
            "the genes are gone;\nwhat remains is $\\bf{how\\ alike}$ the diseases are",
            ha="center", fontsize=9.5, color=GREY)
    ax.set_title("Projected: diseases only", fontsize=11.5, weight="bold")
    ax.set_xlim(0.2, 4.3); ax.set_ylim(0.2, 4.0); ax.axis("off")

    fig.tight_layout()
    fig.savefig(OUT / "bipartite_projection.png", dpi=DPI, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print("  wrote bipartite_projection.png")


def ontology_climb():
    """Inheriting a code from an ancestor when the node itself has none.

    The structure below is the real one - triple-negative breast carcinoma has
    three parents, one per receptor, which is exactly why it is called
    triple-negative. Keep it in step with `s1_helpers.ancestors_of` if the graph
    is rebuilt.
    """
    fig, ax = plt.subplots(figsize=(11.8, 6.6))

    box(ax, (5.0, 5.7), 4.4, 0.5, "triple-negative breast carcinoma",
        TYPE_COLOURS["disease"], 9.5)
    ax.text(7.45, 5.7, "no code", ha="left", va="center", fontsize=8.5,
            color="#B04A4A", style="italic")

    # Three parents at depth 1 - the ontology is a DAG, not a tree.
    receptors = [
        (1.85, "estrogen-receptor\nnegative"),
        (5.0, "progesterone-receptor\nnegative"),
        (8.15, "Her2-receptor\nnegative"),
    ]
    for x, label in receptors:
        box(ax, (x, 4.55), 2.75, 0.66, label, TYPE_COLOURS["disease"], 8.5)
        arrow(ax, (5.0, 5.44), (x, 4.9), "")

    chain = [
        (3.4, "breast carcinoma by gene expression profile", 5.0),
        (2.35, "breast carcinoma", 5.0),
        (1.3, "breast cancer", 5.0),
    ]
    for x, label in receptors:
        arrow(ax, (x, 4.2), (5.0, 3.68), "")
    for i, (y, label, x) in enumerate(chain):
        box(ax, (x, y), 4.4, 0.5, label, TYPE_COLOURS["disease"], 9.5)
        if i < len(chain) - 1:
            arrow(ax, (x, y - 0.26), (x, chain[i + 1][0] + 0.26), "")
        if label != "breast cancer":
            ax.text(7.45, y, "no code", ha="left", va="center", fontsize=8.5,
                    color="#B04A4A", style="italic")

    ax.text(5.15, 5.15, "is_a", ha="left", fontsize=8.5, color=GREY, style="italic")
    ax.text(5.15, 3.9, "is_a", ha="left", fontsize=8.5, color=GREY, style="italic")

    box(ax, (5.0, 0.42), 2.4, 0.5, "ICD10CM:C50", TYPE_COLOURS["icd10"], 10)
    arrow(ax, (5.0, 1.04), (5.0, 0.68), "")
    ax.text(5.15, 0.86, "maps_to", ha="left", va="center", fontsize=8.5,
            color=GREY, style="italic")

    ax.text(-1.5, 4.55, "3 parents,\none per receptor -\nan ontology is a\nnetwork, not a tree",
            ha="left", va="center", fontsize=8.5, color=GREY, style="italic")
    ax.text(5.0, -0.25,
            "4 is_a hops up. The code is real, but it describes\n"
            "the ancestor - not the subtype we started from.",
            ha="center", fontsize=9, color=GREY)

    finish(fig, ax, "ontology_climb.png", (-1.6, 9.6), (-0.75, 6.3),
           "When a disease has no ICD-10 code of its own")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    print(f"Writing diagrams to {OUT}/")
    kg_schema()
    curated_vs_computed()
    bipartite_projection()
    ontology_climb()
    print("done")
