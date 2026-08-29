"""
Session 1 helper functions.

Most of these are adapted from an earlier gene co-expression network workshop,
keeping the same names so they stay recognisable. The originals were written for
co-expression networks - one node type, every edge carrying a
correlation weight, always undirected. A curated knowledge graph breaks all
three of those assumptions, so the ported versions here additionally handle:

  * several node types (gene / disease / icd10) rather than just genes
  * edges that carry no weight at all (`is_a`, `maps_to`)
  * directed graphs, where connected components and clustering differ

Deliberately lightweight: matplotlib, pandas, networkx, numpy and seaborn only
(numpy arrives with pandas anyway). The original module pulled in torch, dgl and
astropy at import time, which is a slow and fragile thing to ask of a workshop
room.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

DATA_DIR = Path(__file__).parent / "../../data/session-1-data/"

# Node types get fixed colours so a gene is the same colour in every figure of
# the session. Anything unrecognised falls through to grey.
TYPE_COLOURS = {
    "gene": "#4C72B0",
    "disease": "#DD8452",
    "icd10": "#55A868",
}
DEFAULT_COLOUR = "#BBBBBB"

MARKER_SHAPES = ["o", "^", "s", "p", "h", "H", "8", "d", "D", "v", "<", ">", "P", "*", "X"]


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------

def load_kg(data_dir: Path | str = DATA_DIR, directed: bool = False) -> nx.Graph:
    """
    Build the Session 1 knowledge graph from the generated CSVs.

    Node attributes: `type`, `name`, `extra`.
    Edge attributes: `type`, and `weight` / `evidence` on `associated_with` only.

    The underlying edges are directional (gene -> disease, disease -> parent,
    disease -> icd10). Undirected is the default because degree, hubs and
    connected components all behave more intuitively that way; pass
    `directed=True` when the direction is the point.
    """
    data_dir = Path(data_dir)
    nodes = pd.read_csv(data_dir / "kg_nodes.csv")
    edges = pd.read_csv(data_dir / "kg_edges.csv")

    G = nx.DiGraph() if directed else nx.Graph()

    for row in nodes.itertuples(index=False):
        G.add_node(
            row.id,
            type=row.type,
            name=row.name if pd.notna(row.name) else row.id,
            extra=row.extra if pd.notna(row.extra) else "",
        )

    for row in edges.itertuples(index=False):
        attrs = {"type": row.type}
        # An undirected graph does not remember which way round an edge was
        # added, so record the parent explicitly. Without this, climbing the
        # ontology silently walks downwards for some edges and upwards for
        # others, depending on networkx's internal adjacency order.
        if row.type == "is_a":
            attrs["parent"] = row.target
        # is_a and maps_to edges are unweighted; leave the attribute off entirely
        # rather than inventing a 0, which would distort every weight plot.
        if pd.notna(row.weight):
            attrs["weight"] = float(row.weight)
        if pd.notna(row.evidence):
            attrs["evidence"] = int(row.evidence)
        G.add_edge(row.source, row.target, **attrs)

    return G


def nodes_of_type(G: nx.Graph, node_type: str) -> list:
    """All node ids of a given `type`."""
    return [n for n, d in G.nodes(data=True) if d.get("type") == node_type]


def edges_of_type(G: nx.Graph, edge_type: str) -> list:
    """All edges of a given `type`, as (u, v) pairs."""
    return [(u, v) for u, v, d in G.edges(data=True) if d.get("type") == edge_type]


def gene_by_symbol(G: nx.Graph, symbol: str) -> str:
    """
    Find a gene node's Ensembl id from its symbol, e.g. "BRCA1" -> ENSG00000012048.

    Participants think in symbols and the graph is keyed by Ensembl id, so this
    saves a manual lookup every time. Raises rather than returning None: a typo
    should stop the cell, not silently produce an empty result further down.
    """
    matches = [n for n, d in G.nodes(data=True)
               if d.get("type") == "gene" and d.get("name") == symbol]
    if not matches:
        raise KeyError(f"No gene node with symbol {symbol!r} in this graph.")
    return matches[0]


# --------------------------------------------------------------------------
# Inferred (co-expression) networks
#
# The counterpart to load_kg. Everything above builds a graph by reading recorded
# facts; everything here builds one by computing correlations from measurements.
# Part 1 puts the two side by side, which is the whole reason this file exists.
# --------------------------------------------------------------------------

def load_expression(data_dir: Path | str = DATA_DIR) -> pd.DataFrame:
    """
    Load the committed expression matrix: patients (rows) x genes (columns).

    This is the TCGA-BRCA transcriptomics view used in Session 2, cut down to the
    genes that are also nodes in the knowledge graph (737 of 760) so that both
    networks describe the same genes. Values are log2-scale and already
    library-size normalised, so correlations can be taken directly.

    Generated by `data-prep/build_coexpression_data.py`.
    """
    path = Path(data_dir) / "coexpr_expression.csv.gz"
    return pd.read_csv(path, index_col=0)


def correlation_matrix(expression: pd.DataFrame) -> pd.DataFrame:
    """
    Gene x gene Pearson correlation, with the diagonal zeroed.

    The diagonal is every gene's perfect correlation with itself. Left in, it
    dominates any "strongest partners" ranking, so it goes to zero here once
    rather than being special-cased at every call site.
    """
    correlations = expression.corr()
    for gene in correlations.index:
        correlations.at[gene, gene] = 0.0
    return correlations


def coexpression_network(expression: pd.DataFrame, threshold: float,
                         absolute: bool = True,
                         correlations: pd.DataFrame | None = None) -> nx.Graph:
    """
    Build an inferred network: an edge wherever two genes correlate above `threshold`.

    Every node is a gene - one node type, unlike the knowledge graph - and every
    edge carries the correlation as its `weight`. Nodes are added for all genes,
    including ones that end up with no edges, because "which genes dropped out at
    this threshold?" is a question worth being able to ask.

    `absolute` thresholds on |r|, so strong anti-correlation counts as a
    relationship. Pass `correlations` to reuse a matrix already computed - it is
    by far the expensive part, and the threshold sweep calls this repeatedly.
    """
    if correlations is None:
        correlations = correlation_matrix(expression)

    genes = list(correlations.index)
    G = nx.Graph()
    for gene in genes:
        G.add_node(gene, type="gene", name=gene)

    values = correlations.to_numpy()
    # Upper triangle only: the matrix is symmetric, so the lower half would just
    # re-add every edge. Done with numpy rather than a double loop because the
    # threshold sweep rebuilds the network once per threshold, and 271,000 pairs
    # in Python is a visible pause in a live session.
    rows, cols = np.triu_indices(len(genes), k=1)
    scores = values[rows, cols]
    keep = (np.abs(scores) if absolute else scores) >= threshold

    G.add_edges_from(
        (genes[i], genes[j], {"type": "coexpressed_with", "weight": float(r)})
        for i, j, r in zip(rows[keep], cols[keep], scores[keep])
    )
    return G


def coexpression_threshold_sweep(expression: pd.DataFrame, thresholds,
                                 correlations: pd.DataFrame | None = None) -> pd.DataFrame:
    """
    Rebuild the network at several thresholds and report what survives each one.

    The point of the table is that there is no principled place to stop. Every row
    is a defensible network built from identical data, and the choice of row is
    the analyst's, not the data's.
    """
    if correlations is None:
        correlations = correlation_matrix(expression)

    rows = []
    for threshold in thresholds:
        G = coexpression_network(expression, threshold, correlations=correlations)
        connected = [n for n, degree in G.degree() if degree > 0]
        components = [c for c in nx.connected_components(G) if len(c) > 1]
        rows.append({
            "threshold": threshold,
            "edges": G.number_of_edges(),
            "connected_genes": len(connected),
            "isolated_genes": G.number_of_nodes() - len(connected),
            "largest_component": max((len(c) for c in components), default=0),
        })
    return pd.DataFrame(rows)


def coexpression_partners(correlations: pd.DataFrame, gene_id: str,
                          G: nx.Graph | None = None, top_n: int = 10) -> pd.DataFrame:
    """
    The genes most strongly co-expressed with one gene, strongest |r| first.

    Pass the knowledge graph as `G` to get gene symbols alongside the Ensembl ids
    - the ranking is unreadable without them.
    """
    if gene_id not in correlations.index:
        raise KeyError(f"{gene_id} has no expression measured in this matrix.")

    ranked = correlations.loc[gene_id].abs().sort_values(ascending=False).head(top_n)
    rows = [{
        "gene_id": partner,
        "symbol": G.nodes[partner]["name"] if G is not None and partner in G else None,
        "r": round(float(correlations.at[gene_id, partner]), 3),
    } for partner in ranked.index]
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Description
# --------------------------------------------------------------------------

def print_graph_info(G: nx.Graph) -> None:
    """
    Print basic information about a graph.

    Unlike the co-expression original this is safe on directed graphs, and it breaks
    the counts down by node and edge type - on a typed graph the totals alone
    hide most of what is interesting.
    """
    directed = G.is_directed()
    print(f"Number of nodes: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")
    print(f"Graph type: {'directed' if directed else 'undirected'}")

    node_types = pd.Series(
        [d.get("type", "untyped") for _, d in G.nodes(data=True)]
    ).value_counts()
    if len(node_types) > 1:
        print("\nNodes by type:")
        for name, count in node_types.items():
            print(f"  {count:>6,}  {name}")

    edge_types = pd.Series(
        [d.get("type", "untyped") for _, _, d in G.edges(data=True)]
    ).value_counts()
    if len(edge_types) > 1:
        print("Edges by type:")
        for name, count in edge_types.items():
            print(f"  {count:>6,}  {name}")

    self_loops = list(nx.selfloop_edges(G))
    print(f"\nSelf-loops: {len(self_loops)}")
    print(f"Graph density: {nx.density(G):.6f}")

    if directed:
        n_components = nx.number_weakly_connected_components(G)
        print(f"Weakly connected components: {n_components}")
    else:
        n_components = nx.number_connected_components(G)
        print(f"Connected components: {n_components}")
        if n_components:
            largest = max(nx.connected_components(G), key=len)
            print(f"Largest component: {len(largest):,} nodes "
                  f"({len(largest) / G.number_of_nodes():.1%} of the graph)")

    print(f"Average clustering coefficient: {nx.average_clustering(G):.4f}")


def get_highest_degree_nodes(G: nx.Graph, top_n: int = 10) -> pd.DataFrame:
    """
    The `top_n` highest-degree nodes.

    Returns a DataFrame rather than the original's list of tuples, so that the
    readable `name` sits next to the opaque id - `ENSG00000012048` means
    little, `BRCA1` means a lot.
    """
    rows = []
    for node, degree in sorted(G.degree(), key=lambda x: x[1], reverse=True)[:top_n]:
        attrs = G.nodes[node]
        rows.append({
            "id": node,
            "name": attrs.get("name", node),
            "type": attrs.get("type", "untyped"),
            "degree": degree,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Cleaning and sparsification
# --------------------------------------------------------------------------

def clean_graph(G: nx.Graph, degree_threshold: int = 1,
                keep_largest_component: bool = True) -> nx.Graph:
    """
    Remove self-loops, isolates and low-degree nodes; optionally keep only the
    largest component. Directed graphs use weakly connected components.
    """
    G = G.copy()
    G.remove_edges_from(nx.selfloop_edges(G))
    G.remove_nodes_from(list(nx.isolates(G)))

    low_degree = [n for n, d in dict(G.degree()).items() if d < degree_threshold]
    G.remove_nodes_from(low_degree)

    if keep_largest_component and G.number_of_nodes():
        components = (nx.weakly_connected_components(G) if G.is_directed()
                      else nx.connected_components(G))
        G = G.subgraph(max(components, key=len)).copy()

    return G


def remove_by_degree(G: nx.Graph, min_degree: int) -> nx.Graph:
    """Drop nodes whose degree is below `min_degree`."""
    G = G.copy()
    G.remove_nodes_from([n for n, d in dict(G.degree()).items() if d < min_degree])
    return G


def threshold_sparsification(G: nx.Graph, threshold: float,
                             keep_unweighted: bool = True) -> nx.Graph:
    """
    Drop edges whose weight is below `threshold`.

    `keep_unweighted` is the knowledge-graph-specific part: `is_a` and
    `maps_to` edges carry no weight, and silently deleting them (which the
    co-expression version would, by treating a missing weight as 0) would tear the
    disease hierarchy out of the graph. Keep them by default.
    """
    out = G.__class__()
    out.add_nodes_from(G.nodes(data=True))
    for u, v, d in G.edges(data=True):
        if "weight" not in d:
            if keep_unweighted:
                out.add_edge(u, v, **d)
        elif d["weight"] >= threshold:
            out.add_edge(u, v, **d)
    return out


def top_percentage_sparsification(G: nx.Graph, top_percentage: float,
                                  keep_unweighted: bool = True) -> nx.Graph:
    """Keep only the top `top_percentage` % of weighted edges."""
    out = G.__class__()
    out.add_nodes_from(G.nodes(data=True))

    weighted, unweighted = [], []
    for u, v, d in G.edges(data=True):
        (weighted if "weight" in d else unweighted).append((u, v, d))

    weighted.sort(key=lambda e: e[2]["weight"], reverse=True)
    keep = max(1, int(len(weighted) * top_percentage / 100)) if weighted else 0
    out.add_edges_from(weighted[:keep])
    if keep_unweighted:
        out.add_edges_from(unweighted)
    return out


# --------------------------------------------------------------------------
# Visualisation
# --------------------------------------------------------------------------

def visualise_graph(G: nx.Graph, title: str = "Knowledge Graph",
                    colour_by: str | None = "type", figsize=(10, 10),
                    seed: int = 0):
    """
    Quick look at a graph, nodes coloured by an attribute.

    The co-expression version drew every node the same blue, which is the right call
    for a co-expression network and the wrong one here - node type is the
    first thing we want to see.
    """
    plt.figure(figsize=figsize)
    pos = nx.spring_layout(G, k=0.1, seed=seed)

    if colour_by:
        colours = [TYPE_COLOURS.get(G.nodes[n].get(colour_by), DEFAULT_COLOUR)
                   for n in G.nodes()]
    else:
        colours = DEFAULT_COLOUR

    nx.draw_networkx_nodes(G, pos, node_size=50, node_color=colours, alpha=0.85)
    nx.draw_networkx_edges(G, pos, width=0.2, alpha=0.5)

    if colour_by:
        present = {G.nodes[n].get(colour_by) for n in G.nodes()}
        plt.legend(handles=[
            mpatches.Patch(color=TYPE_COLOURS.get(v, DEFAULT_COLOUR), label=str(v))
            for v in sorted(present, key=str)
        ], loc="upper left")

    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    plt.show()


def draw_network_with_node_attrs(G, node_attributes=None, communities=None,
                                 title="Network Visualization", color_attr=None,
                                 shape_attr=None, figsize=(20, 10), layout="spring",
                                 cmap_name="tab20", with_labels=False,
                                 node_size=400, seed=0):
    """
    Draw a graph with nodes coloured and/or shaped by attribute.

    Adapted from the co-expression original with two fixes and one signature change:

    * The original grouped nodes for drawing by shape, but built that grouping
      from `shape_attr`. Called with a `color_attr` and no `shape_attr` the
      grouping came out empty and *no nodes were drawn at all* - only edges.
      Colour-only calls now work.
    * Attributes are read from the graph itself by default. The original
      required a separate `node_attributes` dict of the shape
      {attr_name: {node: value}}, which had to be kept in sync by hand. Pass
      one to override to colour by something not on the graph.
    * An empty attribute no longer divides by zero when building the colormap.
    """
    if G is None or G.number_of_nodes() == 0:
        raise ValueError("The graph is empty or not defined.")

    def attr_value(node, attr):
        if node_attributes and attr in node_attributes:
            return node_attributes[attr].get(node)
        return G.nodes[node].get(attr)

    lookup = dict(node_attributes) if node_attributes else {}
    if communities:
        lookup["community"] = {n: i for i, c in enumerate(communities) for n in c}
        node_attributes = lookup
        color_attr = color_attr or "community"

    nodes = list(G.nodes())

    colour_values = sorted({attr_value(n, color_attr) for n in nodes}, key=str) if color_attr else []
    # Node types keep their fixed colour from TYPE_COLOURS so that a gene is the
    # same blue in every figure of the session. Falling back to a colormap here
    # would recolour by position instead, so a disease-only subgraph would come
    # out blue - the colour genes use everywhere else.
    if colour_values and all(v in TYPE_COLOURS for v in colour_values):
        colour_map = {v: TYPE_COLOURS[v] for v in colour_values}
    else:
        cmap = plt.get_cmap(cmap_name)
        colour_map = {
            v: cmap(i / max(len(colour_values) - 1, 1))
            for i, v in enumerate(colour_values)
        }
    node_colours = ([colour_map[attr_value(n, color_attr)] for n in nodes]
                    if color_attr else [DEFAULT_COLOUR] * len(nodes))

    shape_values = sorted({attr_value(n, shape_attr) for n in nodes}, key=str) if shape_attr else []
    shape_map = {v: MARKER_SHAPES[i % len(MARKER_SHAPES)] for i, v in enumerate(shape_values)}

    plt.figure(figsize=figsize)
    layout_fn = getattr(nx, f"{layout}_layout", nx.spring_layout)
    try:
        pos = layout_fn(G, seed=seed)
    except TypeError:  # layouts such as circular_layout take no seed
        pos = layout_fn(G)

    # Group by shape so each marker can be drawn in its own pass. With no
    # shape_attr this is a single group containing every node, which is the
    # bug the original had.
    if shape_attr:
        groups = {s: [n for n in nodes if shape_map[attr_value(n, shape_attr)] == s]
                  for s in shape_map.values()}
    else:
        groups = {"o": nodes}

    index = {n: i for i, n in enumerate(nodes)}
    for shape, group in groups.items():
        if not group:
            continue
        nx.draw_networkx_nodes(
            G, pos,
            nodelist=group,
            node_color=[node_colours[index[n]] for n in group],
            node_shape=shape,
            node_size=node_size,
            edgecolors="white",
            linewidths=0.5,
        )

    nx.draw_networkx_edges(G, pos, width=0.5, alpha=0.5)
    if with_labels:
        nx.draw_networkx_labels(G, pos, font_size=10)

    legends = []
    if shape_attr:
        legends.append(plt.legend(
            handles=[Line2D([0], [0], marker=shape_map[v], color="w", label=str(v),
                            markerfacecolor="k", markersize=10) for v in shape_values],
            title=f"{shape_attr} (shape)", loc="upper left",
            bbox_to_anchor=(1, 0.5), fontsize=12, title_fontsize=12))
    if color_attr:
        legends.append(plt.legend(
            handles=[mpatches.Patch(facecolor=colour_map[v], label=str(v))
                     for v in colour_values],
            title=f"{color_attr} (colour)", loc="upper left",
            bbox_to_anchor=(1, 1), fontsize=12, title_fontsize=12))
    for extra in legends[:-1]:
        plt.gca().add_artist(extra)

    plt.title(title, fontsize=18)
    plt.axis("off")
    plt.tight_layout()
    plt.show()


def gen_graph_legend(G: nx.Graph, attr: str = "type") -> list:
    """
    Legend patches for a graph coloured by `attr`.

    The original took a parallel series of colours and zipped it against
    the attribute values, which relied on the two being in the same order.
    This derives both from the graph, so they cannot disagree.
    """
    values = sorted({d.get(attr) for _, d in G.nodes(data=True)}, key=str)
    return [mpatches.Patch(color=TYPE_COLOURS.get(v, DEFAULT_COLOUR), label=str(v))
            for v in values]


def plot_degree_distribution(G: nx.Graph, bins: int = 30, by_type: bool = False):
    """Degree distribution, optionally split by node type."""
    plt.figure(figsize=(10, 6))
    if by_type:
        frame = pd.DataFrame([
            {"degree": d, "type": G.nodes[n].get("type", "untyped")}
            for n, d in G.degree()
        ])
        sns.histplot(frame, x="degree", hue="type", bins=bins,
                     palette=TYPE_COLOURS, multiple="stack", edgecolor="black")
    else:
        sns.histplot([d for _, d in G.degree()], bins=bins, edgecolor="black")
    plt.title("Degree Distribution")
    plt.xlabel("Degree")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.show()


def visualise_edge_weight_distribution(G: nx.Graph, bins: int = 30):
    """
    Distribution of edge weights.

    Only edges that actually carry a weight are plotted; the original
    indexed `G[u][v]['weight']` directly and would raise a KeyError on the
    first `is_a` edge it met.
    """
    weights = [d["weight"] for _, _, d in G.edges(data=True) if "weight" in d]
    if not weights:
        print("No weighted edges in this graph.")
        return

    unweighted = G.number_of_edges() - len(weights)
    plt.figure(figsize=(10, 6))
    sns.histplot(weights, bins=bins)
    plt.title(f"Distribution of Edge Weights  "
              f"({len(weights):,} weighted, {unweighted:,} unweighted edges)")
    plt.xlabel("Edge Weight")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.show()


# --------------------------------------------------------------------------
# Bipartite projection
# --------------------------------------------------------------------------

def genes_for_disease(G: nx.Graph, disease_id: str) -> set:
    """The set of gene nodes linked to a disease by an `associated_with` edge."""
    if disease_id not in G:
        return set()
    return {n for n in G.neighbors(disease_id) if G.nodes[n].get("type") == "gene"}


def shared_gene_projection(G: nx.Graph, disease_ids=None,
                           min_shared: int = 1) -> nx.Graph:
    """
    Project the bipartite gene-disease graph down onto diseases alone.

    Two diseases are joined if they share at least `min_shared` genes. The edge
    carries both the raw count (`n_shared`) and the Jaccard index (`weight`),
    because the two tell different stories: a well-studied disease shares many
    genes with everything simply by having many genes.

    This is the standard way to turn a two-mode network into a one-mode one,
    and it is where "which diseases resemble each other?" becomes a question
    about graph structure rather than about biology directly.
    """
    if disease_ids is None:
        disease_ids = nodes_of_type(G, "disease")

    gene_sets = {d: genes_for_disease(G, d) for d in disease_ids}
    P = nx.Graph()
    for disease in disease_ids:
        P.add_node(disease, **G.nodes[disease])

    for a, b in itertools.combinations(disease_ids, 2):
        shared = gene_sets[a] & gene_sets[b]
        if len(shared) >= max(min_shared, 1):
            union = gene_sets[a] | gene_sets[b]
            P.add_edge(a, b,
                       n_shared=len(shared),
                       weight=len(shared) / len(union) if union else 0.0,
                       shared_genes=sorted(G.nodes[g]["name"] for g in shared))
    return P


def shared_genes_between(G: nx.Graph, a: str, b: str) -> list:
    """Gene symbols shared by two diseases, for eyeballing a single pair."""
    return sorted(G.nodes[g]["name"] for g in genes_for_disease(G, a) & genes_for_disease(G, b))


# --------------------------------------------------------------------------
# Ontology traversal
# --------------------------------------------------------------------------

def _parent_map(G: nx.Graph) -> dict:
    """child -> {parents}, read from the `parent` stamped on each is_a edge."""
    parents: dict = {}
    for u, v, d in G.edges(data=True):
        if d.get("type") != "is_a":
            continue
        parent = d.get("parent", v)
        child = u if parent == v else v
        parents.setdefault(child, set()).add(parent)
    return parents


def ancestors_of(G: nx.Graph, node: str, with_depth: bool = False) -> list:
    """
    All ancestors reachable by following `is_a` edges upwards, nearest first.

    Breadth-first and alphabetically ordered within each level, so the result is
    identical on every machine and every run. That matters more than it sounds:
    a disease can have several parents, and iterating a set gives a different
    order in each Python process, so an un-sorted version prints a different
    answer every time the notebook is run.

    With `with_depth`, yields (ancestor, depth) where depth is the number of
    `is_a` hops from `node` - a real distance up the tree, not a position in a
    traversal.
    """
    parents = _parent_map(G)
    seen = {node}
    order: list = []
    frontier = [(p, 1) for p in sorted(parents.get(node, ()))]

    while frontier:
        current, depth = frontier.pop(0)
        if current in seen:
            continue
        seen.add(current)
        order.append((current, depth) if with_depth else current)
        frontier.extend((p, depth + 1) for p in sorted(parents.get(current, ())))

    return order


def icd10_for_disease(G: nx.Graph, disease_id: str) -> dict:
    """
    Find an ICD-10 code for a disease, climbing the ontology if necessary.

    Returns a dict with the code, the node it was actually found on, and how
    many `is_a` steps that took. `steps == 0` means the disease carries the
    code itself; anything higher means we inherited it from an ancestor, which
    is a weaker claim and should be reported as such.
    """
    def direct(node):
        return [n for n in G.neighbors(node)
                if G.nodes[n].get("type") == "icd10"] if node in G else []

    own = direct(disease_id)
    if own:
        return {"disease": disease_id, "icd10": own, "found_on": disease_id, "steps": 0}

    # Breadth-first, so the first hit is the *closest* coded ancestor rather than
    # whichever one the traversal happened to reach first.
    for ancestor, depth in ancestors_of(G, disease_id, with_depth=True):
        codes = direct(ancestor)
        if codes:
            return {"disease": disease_id, "icd10": codes,
                    "found_on": ancestor, "steps": depth}

    return {"disease": disease_id, "icd10": [], "found_on": None, "steps": None}


# --------------------------------------------------------------------------
# Bridge to the Session 2 omics data
# --------------------------------------------------------------------------

# The PAM50 subtype labels used in Session 2 are themselves diseases in this
# knowledge graph. That is what lets a patient's molecular subtype be looked up
# as a graph node rather than treated as an opaque string.
PAM50_TO_MONDO = {
    "LumA": "MONDO_0021116",    # luminal A breast carcinoma
    "LumB": "MONDO_0021115",    # luminal B breast carcinoma
    "Basal": "MONDO_0004984",   # basal-like breast carcinoma
    "Her2": "MONDO_0006244",    # HER2 positive breast carcinoma
    "Normal": "MONDO_0006324",  # normal breast-like subtype of breast carcinoma
}


def strip_ensembl_version(ids) -> list:
    """
    Turn versioned Ensembl ids into bare ones: ENSG00000012048.23 -> ENSG00000012048.

    The TCGA matrices carry the version suffix; Open Targets does not. Joining
    the two without this step matches exactly nothing, silently.
    """
    return [str(i).split(".")[0] for i in ids]


def load_omics(path, layers=("transcriptomics",)):
    """
    Load the Session 2 omics pickle, keeping only the layers asked for.

    The full file is ~900 MB, most of it the 200,000-probe methylation matrix,
    so the default keeps transcriptomics only. `meta` (the PAM50 subtype per
    patient) is always returned.

    Returns (dict_of_layers, meta_series).
    """
    import pickle

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"No omics file at {path}. This section is optional - the rest of "
            f"the notebook runs without it."
        )

    with open(path, "rb") as handle:
        blob = pickle.load(handle)

    meta = blob.get("meta")
    kept = {name: blob[name] for name in layers if name in blob}
    missing = [name for name in layers if name not in blob]
    if missing:
        print(f"warning: layers not in the file: {missing}")
    return kept, meta


def map_genes_to_kg(G: nx.Graph, gene_ids) -> pd.DataFrame:
    """
    Look a list of gene ids up in the knowledge graph.

    Accepts versioned or unversioned Ensembl ids. Returns one row per input id
    with whether it was found and, if so, its symbol and degree - so an omics
    gene list becomes an entry point into the graph.
    """
    stripped = strip_ensembl_version(gene_ids)
    rows = []
    for original, gene_id in zip(gene_ids, stripped):
        found = gene_id in G and G.nodes[gene_id].get("type") == "gene"
        rows.append({
            "input_id": original,
            "ensembl_id": gene_id,
            "in_kg": found,
            "symbol": G.nodes[gene_id]["name"] if found else None,
            "degree": G.degree(gene_id) if found else 0,
        })
    return pd.DataFrame(rows)


def diseases_for_genes(G: nx.Graph, gene_ids, top_n: int = 10) -> pd.DataFrame:
    """
    Given a gene list, which diseases does it touch?

    This is the smallest useful knowledge-graph query: one hop out from a set
    of genes, counting where we land. It is the same shape as the queries an
    LLM agent will be asked to plan in Sessions 3 and 4.
    """
    stripped = [g for g in strip_ensembl_version(gene_ids)
                if g in G and G.nodes[g].get("type") == "gene"]

    hits = {}
    for gene in stripped:
        for neighbour in G.neighbors(gene):
            if G.nodes[neighbour].get("type") != "disease":
                continue
            entry = hits.setdefault(neighbour, {"genes": [], "scores": []})
            entry["genes"].append(G.nodes[gene]["name"])
            weight = G.edges[gene, neighbour].get("weight")
            if weight is not None:
                entry["scores"].append(weight)

    rows = [{
        "disease_id": disease,
        "name": G.nodes[disease]["name"],
        "n_genes": len(entry["genes"]),
        "mean_score": (sum(entry["scores"]) / len(entry["scores"])
                       if entry["scores"] else None),
        "genes": ", ".join(sorted(entry["genes"])[:6]),
    } for disease, entry in hits.items()]

    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    return frame.sort_values("n_genes", ascending=False).head(top_n).reset_index(drop=True)


# --------------------------------------------------------------------------
# Evidence types
# --------------------------------------------------------------------------

# What each Open Targets evidence type actually claims. The distinction that
# matters most for teaching is the first two: genetic_association is a claim
# about cause, known_drug is a claim about treatment, and the overall
# association score adds them together as though they were the same thing.
DATATYPE_MEANING = {
    "genetic_association": "inherited variants linked to the disease (cause)",
    "somatic_mutation": "mutations acquired in the tumour (cause)",
    "known_drug": "a drug hitting this target is used for the disease (treatment)",
    "affected_pathway": "the gene sits in a pathway implicated in the disease",
    "rna_expression": "the gene is differentially expressed in the disease",
    "literature": "co-mentioned in publications (text mining)",
    "genetic_literature": "genetic claims extracted from publications",
    "animal_model": "a model organism with this gene disrupted shows the phenotype",
}

# Evidence types that support a claim about *causation* rather than treatment,
# co-mention or correlation. Filtering to these is what separates BRCA2 from
# a tubulin in Part 3.
CAUSAL_DATATYPES = ("genetic_association", "somatic_mutation")


def load_evidence(data_dir: Path | str = DATA_DIR) -> pd.DataFrame:
    """
    The gene-disease edges split by kind of evidence.

    One row per (gene, disease, datatype). The `weight` column is the
    association score contributed by that evidence type alone, so the same
    gene-disease pair appears once per datatype supporting it.
    """
    return pd.read_csv(Path(data_dir) / "kg_evidence.csv")


def evidence_for_pair(evidence: pd.DataFrame, gene_id: str, disease_id: str,
                      G: nx.Graph | None = None) -> pd.DataFrame:
    """Break a single gene-disease edge down into the evidence behind it."""
    rows = evidence[(evidence["source"] == gene_id)
                    & (evidence["target"] == disease_id)].copy()
    rows["means"] = rows["datatype"].map(DATATYPE_MEANING)
    if G is not None and gene_id in G:
        rows.insert(0, "gene", G.nodes[gene_id]["name"])
    return rows.sort_values("weight", ascending=False)


def filter_by_datatype(G: nx.Graph, evidence: pd.DataFrame, datatypes,
                       min_score: float = 0.0) -> nx.Graph:
    """
    Keep only the `associated_with` edges supported by particular evidence types.

    Structural edges (`is_a`, `maps_to`) are always kept - they are not the sort
    of claim evidence types apply to, and dropping them would tear out the
    hierarchy exactly as an over-eager threshold does.

    Pass `CAUSAL_DATATYPES` to keep the edges that assert the gene has something
    to do with causing the disease, and watch how many disappear.
    """
    if isinstance(datatypes, str):
        datatypes = (datatypes,)

    subset = evidence[evidence["datatype"].isin(datatypes)]
    if min_score:
        subset = subset[subset["weight"] >= min_score]
    supported = set(zip(subset["source"], subset["target"]))

    out = G.__class__()
    out.add_nodes_from(G.nodes(data=True))
    for u, v, d in G.edges(data=True):
        if d.get("type") != "associated_with":
            out.add_edge(u, v, **d)
        elif (u, v) in supported or (v, u) in supported:
            out.add_edge(u, v, **d)
    return out


def datatype_summary(evidence: pd.DataFrame) -> pd.DataFrame:
    """How much of each kind of evidence the graph rests on."""
    summary = (
        evidence.groupby("datatype")
        .agg(edges=("weight", "size"), mean_score=("weight", "mean"))
        .sort_values("edges", ascending=False)
        .round(3)
        .reset_index()
    )
    summary["means"] = summary["datatype"].map(DATATYPE_MEANING)
    return summary
