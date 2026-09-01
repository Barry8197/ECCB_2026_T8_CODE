"""
Utilities for *From Multi-Omics to Gene-Disease Discovery: Knowledge Graphs and LLM-Augmented Analysis*
tutorial (ECCB 2026, Geneva).

This is the trimmed, **Session 4 mini-challenge only** version of the
tutorial's helper module. The full workshop module also contains Session 1
material (co-expression networks, ontology traversal, evidence-type
filtering) and MOFA-fitting/diagnostic-plotting helpers used while *building*
the pretrained model — none of that is needed to run the mini challenge, so
it has been left out here to keep this file short and easy to read during
the session. If you need those extras for other parts of the workshop, keep
using the full `s4_helpers.py`.

Exported functions
------------------
- load_kg
- print_graph_info
- map_genes_to_kg
- diseases_for_genes
- load_omics
- project_test_patients_to_mofa_factors
- evaluate_predictions
"""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from IPython.display import display

from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    ConfusionMatrixDisplay,
)

__all__ = [
    "load_kg",
    "print_graph_info",
    "map_genes_to_kg",
    "diseases_for_genes",
    "load_omics",
    "project_test_patients_to_mofa_factors",
    "evaluate_predictions",
]

# Default location of the knowledge-graph CSVs (kg_nodes.csv / kg_edges.csv).
DATA_DIR = Path("/data/session-1-data/")


# --------------------------------------------------------------------------
# Knowledge graph: loading and inspection
# --------------------------------------------------------------------------

def load_kg(data_dir: Path | str = DATA_DIR, directed: bool = False) -> nx.Graph:
    """
    Build the knowledge graph from the generated CSVs.

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
        # added, so record the parent explicitly.
        if row.type == "is_a":
            attrs["parent"] = row.target
        # is_a and maps_to edges are unweighted; leave the attribute off
        # entirely rather than inventing a 0.
        if pd.notna(row.weight):
            attrs["weight"] = float(row.weight)
        if pd.notna(row.evidence):
            attrs["evidence"] = int(row.evidence)
        G.add_edge(row.source, row.target, **attrs)

    return G


def print_graph_info(G: nx.Graph) -> None:
    """
    Print basic information about a graph: node/edge counts by type, density,
    self-loops, and connected components. Safe on both directed and
    undirected graphs.
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


# --------------------------------------------------------------------------
# Bridge from an omics gene list into the knowledge graph
# --------------------------------------------------------------------------

def strip_ensembl_version(ids) -> list:
    """
    Turn versioned Ensembl ids into bare ones: ENSG00000012048.23 -> ENSG00000012048.

    The TCGA matrices carry the version suffix; the knowledge graph does not.
    Joining the two without this step matches exactly nothing, silently.
    """
    return [str(i).split(".")[0] for i in ids]


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
    of genes, counting where we land.
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
# Omics loading
# --------------------------------------------------------------------------

def load_omics(
    data_dir: str | Path,
    omic_filename: str,
    omic_keys: Sequence[str],
) -> tuple[dict[str, pd.DataFrame], pd.Series]:
    """
    Load TCGA-BRCA multi-omics views and subtype labels from a pickled bundle.

    The function reads ``{omic_filename}.pkl`` from `data_dir`. The pickle is
    expected to contain one DataFrame per requested omics view key, plus a
    label vector under the key ``"meta"``.

    The function also checks that all requested views share an identical
    patient index.

    Parameters
    ----------
    data_dir : str | pathlib.Path
        Directory containing ``{omic_filename}.pkl``.
    omic_filename : str
        Base filename (without ``.pkl``) of the pickled omics bundle.
    omic_keys : Sequence[str]
        Names of the omics views to load, e.g.
        ``["transcriptomics", "proteomics", "methylation"]``.

    Returns
    -------
    X_views : dict[str, pandas.DataFrame]
        Mapping from each requested view name to a copy of its feature matrix.
    y : pandas.Series
        Copy of the labels (subtype), indexed by patient/sample ID.

    Raises
    ------
    ValueError
        If `omic_keys` is empty, or if patient indices across views are not identical.
    FileNotFoundError
        If ``{omic_filename}.pkl`` does not exist in `data_dir`.
    KeyError
        If any requested omics key (or ``"meta"``) is missing from the pickle.

    Notes
    -----
    Prints basic dataset diagnostics (view dimensions and label counts).
    """
    data_dir = Path(data_dir)

    if omic_keys is None or len(omic_keys) == 0:
        raise ValueError("omic_keys is required and must contain at least one key.")

    # ---- Load the bundle -------------------------------------------------
    omics_path = data_dir / f"{omic_filename}.pkl"
    if not omics_path.exists():
        raise FileNotFoundError(f"{omic_filename}.pkl not found in: {data_dir}")

    omics = pd.read_pickle(omics_path)

    # ---- Validate expected keys -----------------------------------------
    required = list(omic_keys) + ["meta"]
    missing = [k for k in required if k not in omics]
    if missing:
        raise KeyError(f"Missing keys in {omic_filename}.pkl: {missing}")

    # ---- Print dimensions (quick sanity check) --------------------------
    print("Omic view dimensions:")
    for key in omic_keys:
        n_patients, n_features = omics[key].shape
        print(f"  {key:15s}: {n_patients:4d} patients x {n_features:6d} features")

    # ---- Assert alignment across requested views ------------------------
    reference_index = omics[omic_keys[0]].index
    for key in omic_keys[1:]:
        if not reference_index.equals(omics[key].index):
            raise ValueError(
                f"Patient index of '{key}' does not match '{omic_keys[0]}'. "
                "Views must be pre-aligned."
            )

    # ---- Build outputs ---------------------------------------------------
    X_views = {key: omics[key].copy() for key in omic_keys}
    y = omics["meta"].copy()

    # ---- Optional: label distribution (nice in notebooks) ---------------
    print("\nSubtype counts:")
    display(y.value_counts())

    return X_views, y


# --------------------------------------------------------------------------
# MOFA: projecting held-out patients into a fixed factor space
# --------------------------------------------------------------------------

def project_test_patients_to_mofa_factors(model, X_train_by_view, X_test_by_view, train_factors, view_names):
    """
    Project held-out samples into a trained MOFA factor space.

    MOFA is fitted on training samples only. For held-out samples, this
    function:
    - fixes learned weights (W) per view
    - computes a pseudo-inverse-based projection from scaled test features to factors
    - calibrates the raw projection to the trained factor scale using a linear
      map fitted on training samples only
    - averages projected factor values across views

    Parameters
    ----------
    model : Any
        Fitted MOFA model object providing ``get_weights(views=..., df=True)``.
    X_train_by_view : dict[str, pandas.DataFrame]
        Training feature matrices per view, indexed by sample ID.
    X_test_by_view : dict[str, pandas.DataFrame]
        Test feature matrices per view, indexed by sample ID.
    train_factors : pandas.DataFrame
        Training factor matrix (Z) indexed by sample ID with columns being factor names.
    view_names : list[str]
        Ordered list of view names to project with.

    Returns
    -------
    pandas.DataFrame
        Projected test factor values, indexed by test sample ID with the same
        factor columns as `train_factors`.

    Notes
    -----
    - Scaling is done per view using training-set mean and standard deviation.
    - Features are intersected across MOFA weights, train data, and test data
      for safety.
    - Test labels are never used anywhere in this function.
    """
    factor_columns = train_factors.columns.astype(str).tolist()
    projected_test_by_view = []

    for view_name in view_names:
        weights = model.get_weights(views=view_name, df=True)
        weights.columns = weights.columns.astype(str)
        weights = weights.reindex(columns=factor_columns)

        common_features = weights.index.intersection(X_train_by_view[view_name].columns)
        common_features = common_features.intersection(X_test_by_view[view_name].columns)
        weights = weights.loc[common_features]

        X_train_view = X_train_by_view[view_name].loc[:, common_features].astype(float)
        X_test_view = X_test_by_view[view_name].loc[:, common_features].astype(float)
        X_train_view.index = X_train_view.index.astype(str)
        X_test_view.index = X_test_view.index.astype(str)

        train_mean = X_train_view.mean(axis=0)
        train_std = X_train_view.std(axis=0, ddof=0).replace(0, 1)
        X_train_scaled = (X_train_view - train_mean) / train_std
        X_test_scaled = (X_test_view - train_mean) / train_std

        raw_train_projection = X_train_scaled.to_numpy() @ np.linalg.pinv(weights.to_numpy()).T
        raw_test_projection = X_test_scaled.to_numpy() @ np.linalg.pinv(weights.to_numpy()).T

        train_design = np.column_stack([raw_train_projection, np.ones(raw_train_projection.shape[0])])
        test_design = np.column_stack([raw_test_projection, np.ones(raw_test_projection.shape[0])])
        train_target = train_factors.loc[X_train_view.index, factor_columns].to_numpy()
        calibration = np.linalg.lstsq(train_design, train_target, rcond=None)[0]
        projected_values = test_design @ calibration

        projected = pd.DataFrame(projected_values, index=X_test_view.index, columns=factor_columns)
        projected_test_by_view.append(projected)

    return sum(projected_test_by_view) / len(projected_test_by_view)


# --------------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------------

def plot_confusion_matrix(y_test, y_pred) -> None:
    """Plot a confusion matrix for held-out subtype predictions."""
    fig, ax = plt.subplots(figsize=(5, 5))
    ConfusionMatrixDisplay.from_predictions(y_test, y_pred, xticks_rotation=45, ax=ax)
    ax.set_title("Confusion matrix")
    fig.tight_layout()
    plt.show()


def evaluate_predictions(y_true: np.ndarray, y_pred: np.ndarray, title: str) -> dict:
    """
    Compute and display common classification metrics.

    Prints a title header, accuracy, and balanced accuracy, then shows a
    confusion matrix (via :func:`plot_confusion_matrix`).

    Parameters
    ----------
    y_true : numpy.ndarray
        Ground-truth labels of shape `(n_samples,)`.
    y_pred : numpy.ndarray
        Predicted labels of shape `(n_samples,)`.
    title : str
        Title used in printed output and metric labelling.

    Returns
    -------
    dict
        Summary metrics with keys: ``model``, ``accuracy``, ``balanced_accuracy``.
    """
    sep = "-" * len(title)
    accuracy = accuracy_score(y_true, y_pred)
    balanced_accuracy = balanced_accuracy_score(y_true, y_pred)
    print(f"\n{title}\n{sep}")
    print(f"  Accuracy          : {accuracy:.3f} ")
    print(f"  Balanced accuracy : {balanced_accuracy:.3f}")
    print()
    plot_confusion_matrix(y_true, y_pred)

    return {
        "model": title,
        "accuracy": accuracy,
        "balanced_accuracy": balanced_accuracy,
    }