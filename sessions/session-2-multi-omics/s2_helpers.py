"""
Utilities for *From Multi-Omics to Gene–Disease Discovery: Knowledge Graphs and LLM-Augmented Analysis*
tutorial (ECCB 2026, Geneva).

This module contains utility functions used in Session 2 (multi-omics / MOFA)
and a small set of general evaluation helpers.

Exported functions
------------------
- load_omics
- evaluate_predictions
- plot_confusion_matrix
- generate_diagnostic_plots
- build_mofa_matrix_input
- fit_mofa
- select_active_factors
- project_test_patients_to_mofa_factors
- eta_squared_by_factor
- fit_factor_classifier
- plot_r2_heatmap
- plot_factor_boxplots_by_subtype
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
    ConfusionMatrixDisplay
)
from sklearn.linear_model import LogisticRegression

from mofapy2.run.entry_point import entry_point
import mofax as mfx

__all__ = [
    "load_omics",
    "evaluate_predictions",
    "plot_confusion_matrix",
    "generate_diagnostic_plots",
    "build_mofa_matrix_input",
    "fit_mofa",
    "select_active_factors",
    "project_test_patients_to_mofa_factors",
    "eta_squared_by_factor",
    "fit_factor_classifier",
    "plot_r2_heatmap",
    "plot_factor_boxplots_by_subtype",
]

# ---------------------------------------------------------------------------
# SESSION - 2 - LLM Functions
# ---------------------------------------------------------------------------

def build_mofa_matrix_input(X_by_omic):
    """
    Build MOFA input matrices from aligned omics tables.

    `mofapy2` expects a nested structure where each *view* (omics modality) contains
    one or more *groups*, and each entry is a numeric matrix of shape
    `(n_samples, n_features)`.

    This helper assumes:
    - all input DataFrames are already aligned on an identical patient index
    - a single group is used: ``"TCGA-BRCA_train"``

    Parameters
    ----------
    X_by_omic : dict[str, pandas.DataFrame]
        Mapping from view name (e.g. ``"transcriptomics"``) to a 2D feature table
        indexed by patient/sample ID. Each DataFrame has shape
        `(n_samples, n_features_view)`.

    Returns
    -------
    data : list[list[numpy.ndarray]]
        MOFA data structure: ``data[view][group] = X`` where each ``X`` is a
        ``float32`` NumPy array of shape `(n_samples, n_features_view)`.
    view_names : list[str]
        Names of the omics views in the same order as `data`.
    feature_names : list[list[str]]
        Feature names per view in the same order as `view_names`.
    sample_names : list[list[str]]
        Sample IDs per group. For a single group this is `[sample_id_list]`.
    group_names : list[str]
        Group names. For this tutorial this is `["TCGA-BRCA_train"]`.

    Raises
    ------
    AssertionError
        If any view's sample index does not match the first view's sample index.
    """
    data, view_names, feature_names = [], [], []
    mofa_sample_ids = next(iter(X_by_omic.values())).index.astype(str)

    for view_name, X in X_by_omic.items():
        assert X.index.astype(str).equals(mofa_sample_ids), f"Patient index differs in {view_name}"
        view_names.append(view_name)
        feature_names.append(X.columns.astype(str).tolist())
        data.append([X.to_numpy(dtype=np.float32)])  # one group, many views

    sample_names = [mofa_sample_ids.tolist()]
    group_names = ["TCGA-BRCA_train"]

    return data, view_names, feature_names, sample_names, group_names

def fit_mofa(data, view_names, feature_names, sample_names, group_names,
             max_factors, iterations, random_state, outfile):
    """
    Fit a MOFA model (mofapy2) and return the trained model object.

    This function configures a compact MOFA model for continuous (Gaussian)
    multi-omics data, enabling ARD (automatic relevance determination) at both
    the factor and weight level to encourage sparsity/shrinkage of uninformative
    factors.

    Parameters
    ----------
    data : list[list[numpy.ndarray]]
        Nested data structure as produced by :func:`build_mofa_matrix_input`.
    view_names : list[str]
        View names corresponding to entries in `data`.
    feature_names : list[list[str]]
        Feature names per view.
    sample_names : list[list[str]]
        Sample names per group.
    group_names : list[str]
        Group names.
    max_factors : int
        Maximum number of latent factors to learn.
    iterations : int
        Number of training iterations.
    random_state : int
        Random seed passed to MOFA training.
    outfile : str | pathlib.Path
        Path for MOFA's training output and/or saved model.

    Returns
    -------
    model : mofapy2.run.entry_point.entry_point
        Trained MOFA entry point object.

    Notes
    -----
    The function calls ``model.run()`` and then saves the model (if needed).
    """
    model = entry_point()

    model.set_data_options(
        scale_views=True,
        scale_groups=False,
        center_groups=True,
        use_float32=True,
    )

    model.set_data_matrix(
        data=data,
        likelihoods=["gaussian"] * len(view_names),
        views_names=view_names,
        groups_names=group_names,
        samples_names=sample_names,
        features_names=feature_names,
    )

    model.set_model_options(
        factors=max_factors,
        spikeslab_factors=False,
        spikeslab_weights=True,
        ard_factors=True,
        ard_weights=True,
    )

    model.set_train_options(
        iter=iterations,
        convergence_mode="fast",
        seed=random_state,
        verbose=False,
        quiet=True,
        outfile=str(outfile),
    )

    model.build()
    model.run()

    if not Path(outfile).exists():
        model.save(outfile=str(outfile))

    return model

def select_active_factors(mofa_model_mfx, min_total_r2, max_factors):
    """
    Select "active" MOFA factors using variance explained (R2).

    This function retrieves the variance-explained table from a fitted MOFA model,
    sums each factor's R2 across all views, and keeps factors whose total R2 meets
    a minimum threshold. If none pass, it falls back to the single best factor.

    Parameters
    ----------
    mofa_model_mfx : mofax.core.mofa_model.MofaModel | Any
        A `mofax`-wrapped MOFA model providing ``get_r2()``.
    min_total_r2 : float
        Minimum total (summed across views) R2 required to keep a factor.
    max_factors : int
        Maximum number of factors considered/trained (kept for API compatibility /
        notebook parity).

    Returns
    -------
    active_factor_cols : list[str]
        Selected factor names (e.g. ``["Factor1", "Factor3"]``).
    factor_r2_summary : pandas.DataFrame
        Table with columns:
        - ``factor``: factor name
        - ``total_r2``: summed R2 across views
        Sorted descending by ``total_r2``.

    Notes
    -----
    `max_factors` is not used directly in the current implementation but is kept
    to match the tutorial’s calling signature.
    """
    r2_all = mofa_model_mfx.get_r2().rename(
        columns={"Factor": "factor", "View": "view", "Group": "group_mofax", "R2": "r2"}
    )

    factor_r2_summary = (
        r2_all.groupby("factor", as_index=False)["r2"].sum()
        .rename(columns={"r2": "total_r2"})
        .sort_values("total_r2", ascending=False)
    )

    active_factor_summary = factor_r2_summary[factor_r2_summary["total_r2"] >= min_total_r2].copy()
    if active_factor_summary.empty:
        active_factor_summary = factor_r2_summary.head(1).copy()

    active_factor_cols = active_factor_summary["factor"].tolist()
    return active_factor_cols, factor_r2_summary

def project_test_patients_to_mofa_factors(model, X_train_by_view, X_test_by_view, train_factors, view_names):
    """
    Project held-out samples into a trained MOFA factor space. 
    
    MOFA is fitted on training samples only. For held-out samples, this function:
    - fixes learned weights (W) per view
    - computes a pseudo-inverse-based projection from scaled test features to factors
    - calibrates the raw projection to the trained factor scale using a linear map fitted on training samples only
    - averages projected factor values across views

    Parameter
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
    - Features are intersected across MOFA weights, train data, and test data for safety.
    
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

def eta_squared_by_factor(factor_table, labels):
    """
    Compute one-way ANOVA eta-squared ($ \\\\eta^2 $) for each factor.

    Eta-squared here is interpreted as the fraction of *factor variance* explained
    by group (subtype) membership, computed without fitting a predictive model.

    Parameters
    ----------
    factor_table : pandas.DataFrame
        Factor values (e.g. MOFA Z) with shape `(n_samples, n_factors)`,
        indexed by sample ID.
    labels : pandas.Series
        Group labels (e.g. subtype) indexed by the same sample IDs as `factor_table`.

    Returns
    -------
    pandas.DataFrame
        Sorted table with columns:
        - ``factor``: factor column name
        - ``eta_squared``: eta-squared value (float)

    Notes
    -----
    - If a factor has zero total variance, eta-squared is reported as ``NaN``.
    - Labels are coerced to string to avoid mixed label types.
    """
    rows = []
    labels = labels.astype(str)

    for factor in factor_table.columns:
        values = factor_table[factor]
        grand_mean = values.mean()
        ss_total = ((values - grand_mean) ** 2).sum()
        ss_between = 0.0

        for _, idx in labels.groupby(labels).groups.items():
            group_values = values.loc[idx]
            ss_between += len(group_values) * (group_values.mean() - grand_mean) ** 2

        eta2 = ss_between / ss_total if ss_total > 0 else np.nan
        rows.append({"factor": factor, "eta_squared": eta2})

    return pd.DataFrame(rows).sort_values("eta_squared", ascending=False)

def fit_factor_classifier(factors_df, y, train_ids, test_ids, model_name):
    """
    Train a logistic regression classifier on factor values and evaluate on a held-out set.

    This is a lightweight diagnostic answering: *do the learned latent factors
    preserve subtype information?*

    Parameters
    ----------
    factors_df : pandas.DataFrame
        Factor matrix (Z) indexed by sample ID, columns are factor names.
    y : pandas.Series
        True labels indexed by sample ID.
    train_ids : array-like
        Sample IDs to use for training.
    test_ids : array-like
        Sample IDs to use for evaluation.
    model_name : str
        Label used in metric output.

    Returns
    -------
    clf : sklearn.linear_model.LogisticRegression
        Fitted classifier.
    pred : numpy.ndarray
        Predicted labels for `test_ids`.
    metrics : dict
        Dictionary returned by :func:`evaluate_predictions`.
    """
    clf = LogisticRegression(
        max_iter=2000,
        class_weight="balanced",
        random_state=42,
    )
    clf.fit(factors_df.loc[train_ids], y.loc[train_ids])
    pred = clf.predict(factors_df.loc[test_ids])
    metrics = evaluate_predictions(y.loc[test_ids], pred, model_name)
    return clf, pred, metrics

def plot_r2_heatmap(r2_all, view_names, active_factor_cols, output_path):
    """
    Plot an annotated heatmap of MOFA variance explained (R2) and save to disk.

    Parameters
    ----------
    r2_all : pandas.DataFrame
        Long-form R2 table with columns:
        - ``factor``
        - ``view``
        - ``r2`` (numeric)
        Typically derived from ``mofa_model_mfx.get_r2()`` with column renaming.
    view_names : list[str]
        View order for plotting (rows).
    active_factor_cols : list[str]
        Factor order for plotting (columns).
    output_path : str | pathlib.Path
        Output path for a PNG file.

    Returns
    -------
    None
        Figure is saved and closed.
    """
    r2_heatmap = (
        r2_all[r2_all["factor"].isin(active_factor_cols)]
        .pivot(index="view", columns="factor", values="r2")
        .reindex(index=view_names, columns=active_factor_cols)
    )

    fig, ax = plt.subplots(figsize=(1.1 * len(active_factor_cols) + 3, 3.6))
    im = ax.imshow(r2_heatmap, aspect="auto", cmap="Blues")

    ax.set_title("MOFA R2 heatmap: views x active factors")
    ax.set_xlabel("MOFA factor")
    ax.set_ylabel("Omics view")
    ax.set_xticks(range(r2_heatmap.shape[1]))
    ax.set_xticklabels(r2_heatmap.columns, rotation=45, ha="right")
    ax.set_yticks(range(r2_heatmap.shape[0]))
    ax.set_yticklabels(r2_heatmap.index)

    for row in range(r2_heatmap.shape[0]):
        for col in range(r2_heatmap.shape[1]):
            value = r2_heatmap.iloc[row, col]
            if pd.notna(value):
                ax.text(col, row, f"{value:.1f}", ha="center", va="center", fontsize=8)

    fig.colorbar(im, ax=ax, label="MOFA R2 (%)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

def plot_factor_boxplots_by_subtype(factors_df, train_ids, y_train, top_factors, output_path):
    """
    Plot boxplots of selected factors stratified by subtype and save to disk.

    Parameters
    ----------
    factors_df : pandas.DataFrame
        Factor matrix indexed by sample ID.
    train_ids : array-like
        Training sample IDs to include in the plot.
    y_train : pandas.Series
        Training labels (subtypes).
    top_factors : list[str]
        Factor column names to plot (one panel per factor).
    output_path : str | pathlib.Path
        Output path for a PNG file.

    Returns
    -------
    None
        Figure is saved and closed.
    """
    train_factors_for_plot = factors_df.loc[train_ids.astype(str), top_factors]
    train_labels_for_plot = y_train.astype(str)
    subtypes = train_labels_for_plot.dropna().unique()

    n_panels = len(top_factors)
    fig, axes = plt.subplots(1, n_panels, figsize=(3.75 * n_panels, 4), squeeze=False)

    for ax, factor in zip(axes.ravel(), top_factors):
        groups = [
            train_factors_for_plot.loc[train_labels_for_plot == subtype, factor].dropna()
            for subtype in subtypes
        ]
        ax.boxplot(groups, tick_labels=subtypes, showfliers=False)
        ax.set_title(factor)
        ax.set_ylabel("Factor value")
        ax.tick_params(axis="x", rotation=45)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

def plot_ranked_feature_weights(mofa_model_mfx, ranked_factors, view, output_path, n_features=5):
    """
    Plot ranked, signed MOFA feature weights for a given view and save to disk.

    For each requested factor, this calls :func:`mofax.plot_weights_ranked` and
    saves a multi-panel figure.

    Parameters
    ----------
    mofa_model_mfx : mofax.core.mofa_model.MofaModel | Any
        `mofax` model wrapper providing plotting helpers.
    ranked_factors : list[str]
        Factor names to plot.
    view : str
        View name to plot weights for.
    output_path : str | pathlib.Path
        Output path for a PNG file.
    n_features : int, default=5
        Number of features to label at each extreme of the ranked plot.

    Returns
    -------
    None
        Figure is saved and closed.
    """
    fig, axes = plt.subplots(1, len(ranked_factors), figsize=(5 * len(ranked_factors), 4.5), squeeze=False)

    for ax, ranked_factor in zip(axes.ravel(), ranked_factors):
        mfx.plot_weights_ranked(
            mofa_model_mfx,
            factor=ranked_factor,
            view=view,
            n_features=n_features,
            ax=ax,
        )
        ax.set_title(f"{view}: {ranked_factor}")

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

def generate_diagnostic_plots(mofa_model_mfx, factors_df, factor_r2_summary, view_names, active_factor_cols, factor_subtype_assoc, train_ids, y_train, y_test, mofa_pred, output_dir, top_view_for_weights="transcriptomics"):
    """
    Generate and save the core MOFA diagnostic plots used in the tutorial notebook.

    This is a single entry point that produces several complementary diagnostics:

    1. **R2 heatmap**: MOFA variance explained per (view, factor).
    2. **Factor boxplots**: distributions of top subtype-associated factors by subtype
       (training patients only).
    3. **Confusion matrix**: held-out subtype prediction errors from a classifier
       trained on factor values.
    4. **Ranked feature weights**: strongest positive/negative feature loadings in a
       chosen view for the most subtype-associated factors.

    Parameters
    ----------
    mofa_model_mfx : mofax.core.mofa_model.MofaModel | Any
        `mofax` model wrapper providing ``get_r2()`` and plotting.
    factors_df : pandas.DataFrame
        Factor matrix (Z) for all samples, indexed by sample ID.
    factor_r2_summary : pandas.DataFrame
        Factor-level R2 summary table (kept for notebook parity; not directly used).
    view_names : list[str]
        View names (row order for R2 heatmap).
    active_factor_cols : list[str]
        Selected active factor names (column order for R2 heatmap).
    factor_subtype_assoc : pandas.DataFrame
        Table ranking factors by subtype association (e.g. output of
        :func:`eta_squared_by_factor`) with a ``factor`` column.
    train_ids : array-like
        Training sample IDs.
    y_train : pandas.Series
        Training labels.
    y_test : pandas.Series
        Test labels.
    mofa_pred : array-like
        Predicted labels for the test set (for confusion matrix).
    output_dir : str | pathlib.Path
        Directory where PNGs will be written.
    top_view_for_weights : str, default="transcriptomics"
        View to use when plotting ranked feature weights.

    Returns
    -------
    None
        Side effect: saves PNG files under `output_dir`.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Recreate the "all views x all factors" R2 table used by the heatmap.
    r2_all = mofa_model_mfx.get_r2().rename(
        columns={"Factor": "factor", "View": "view", "Group": "group_mofax", "R2": "r2"}
    )

    plot_r2_heatmap(
        r2_all, view_names, active_factor_cols,
        output_dir / "part2_mofa_r2_heatmap.png",
    )

    top_factors = factor_subtype_assoc.head(4)["factor"].tolist()

    plot_factor_boxplots_by_subtype(
        factors_df, train_ids, y_train, top_factors,
        output_dir / "part2_mofa_factor_boxplots.png",
    )

    plot_confusion_matrix(
        y_test, mofa_pred
    )

    plot_ranked_feature_weights(
        mofa_model_mfx, top_factors[:3], top_view_for_weights,
        output_dir / "part2_mofa_ranked_weights.png",
    )

    print(f"Saved diagnostic plots to: {output_dir}")

# ---------------------------------------------------------------------------
# SESSION - 2 - General Functions
# ---------------------------------------------------------------------------

def plot_confusion_matrix(y_test, y_pred):
    """
    Plot a confusion matrix for held-out subtype predictions.

    Parameters
    ----------
    y_test : array-like
        True labels for the evaluation set.
    y_pred : array-like
        Predicted labels for the evaluation set.

    Returns
    -------
    None
        Displays the plot via Matplotlib.
    """
    fig, ax = plt.subplots(figsize=(5, 5))
    ConfusionMatrixDisplay.from_predictions(y_test, y_pred, xticks_rotation=45, ax=ax)
    ax.set_title("Confusion matrix")
    fig.tight_layout()
    plt.show()

def evaluate_predictions(y_true: np.ndarray, y_pred: np.ndarray, title: str) -> None:
    """
    Compute and display common classification metrics.

    The function prints:
    - a title header
    - accuracy
    - balanced accuracy
    - a confusion matrix (via :func:`plot_confusion_matrix`)

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
        Summary metrics with keys:
        - ``model``
        - ``accuracy``
        - ``balanced_accuracy``

    Notes
    -----
    `classification_report` is imported but not printed in the current implementation.
    """
    sep = "─" * len(title)
    accuracy = accuracy_score(y_true, y_pred)
    balanced_accuracy =  balanced_accuracy_score(y_true, y_pred)
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

def load_omics(
    data_dir: str | Path,
    omic_keys: Sequence[str],
    ) -> tuple[dict[str, pd.DataFrame], pd.Series]:
    """
    Load TCGA-BRCA multi-omics views and subtype labels from a pickled bundle.

    The function reads ``omics.pkl`` from `data_dir`. The pickle is expected to contain:
    - one DataFrame per requested omics view key
    - a label vector under the key ``"meta"``

    The function also checks that all requested views share an identical patient index.

    Parameters
    ----------
    data_dir : str | pathlib.Path
        Directory containing ``omics.pkl``.
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
        If `omic_keys` is empty.
    FileNotFoundError
        If ``omics.pkl`` does not exist in `data_dir`.
    KeyError
        If any requested omics key (or ``"meta"``) is missing from the pickle.
    ValueError
        If patient indices across views are not identical.

    Notes
    -----
    The function prints basic dataset diagnostics (view dimensions and label counts).
    """
    data_dir = Path(data_dir)

    if omic_keys is None or len(omic_keys) == 0:
        raise ValueError("omic_keys is required and must contain at least one key.")

    # ---- Load the bundle -------------------------------------------------
    omics_path = data_dir / "omics.pkl"
    if not omics_path.exists():
        raise FileNotFoundError(f"omics.pkl not found in: {data_dir}")

    omics = pd.read_pickle(omics_path)

    # ---- Validate expected keys -----------------------------------------
    required = list(omic_keys) + ["meta"]
    missing = [k for k in required if k not in omics]
    if missing:
        raise KeyError(f"Missing keys in omics.pkl: {missing}")

    # ---- Print dimensions (quick sanity check) --------------------------
    print("Omic view dimensions:")
    for key in omic_keys:
        n_patients, n_features = omics[key].shape
        print(f"  {key:15s}: {n_patients:4d} patients × {n_features:6d} features")

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
