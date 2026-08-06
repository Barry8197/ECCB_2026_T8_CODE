"""
Part 2 -- Correlation-Based / Factor-Based Multi-Omic Integration with MOFA
Minimal working example (.py version)

This script is a condensed version of the original notebook. It keeps the core
computational steps plus a compact set of the notebook's diagnostic plots:

1. Load pre-aligned omics tables (transcriptomics, proteomics, methylation) + subtype labels.
2. Create one shared patient-level train/test split and select variable features
   for high-dimensional views.
3. Build the nested MOFA input structure and fit a MOFA model with mofapy2.
4. Read training factor values and project held-out test patients into the same
   factor space (no refitting).
5. Select "active" factors from MOFA's variance-explained (R2) table.
6. Quantify factor <-> subtype association with eta-squared.
7. Train a simple logistic regression classifier on MOFA factors and evaluate it
   on the held-out test patients.
8. Save metrics, predictions, and factor tables to outputs/.
9. Generate and save a compact set of diagnostic plots (R2 heatmap, factor
   scatter/boxplots by subtype, confusion matrix, ranked feature weights) via
   a single `generate_diagnostic_plots` wrapper function.

Many discussion cells and the full set of exploratory `mofax` plotting
functions from the original notebook have been trimmed to keep this a
minimal, runnable example; `generate_diagnostic_plots` covers the highest-
value plot from each section instead of every variant shown in the notebook.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    ConfusionMatrixDisplay,
)
from sklearn.model_selection import train_test_split

from mofapy2.run.entry_point import entry_point
import mofax as mfx


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
RANDOM_STATE = 42
DATA_DIR = Path("../../data_tmp/TCGA-BRCA")
OUTPUT_DIR = Path("outputs").resolve()

TEST_SIZE = 0.25
N_TOP_VARIABLE_HIGH_DIM_FEATURES = 2000
HIGH_DIMENSIONAL_VIEWS = ["transcriptomics", "methylation", "proteomics"]
MAX_FACTORS = 10
MIN_TOTAL_R2 = 2.5
MOFA_ITERATIONS = 500
SAVE_MOFA_HDF5 = True
LOAD_MOFA_HDF5 = True
MOFA_HDF5_FILE = OUTPUT_DIR / (
    f"trained_mofaplus_train_var{N_TOP_VARIABLE_HIGH_DIM_FEATURES}"
    f"_max{MAX_FACTORS}_ard_model.hdf5"
)


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------
def load_omics_data(data_dir):
    """Load the pre-aligned multi-omics pickle and split it into X (views) and y (labels).

    Expects a dict with keys 'transcriptomics', 'proteomics', 'methylation', and
    'meta' (the subtype labels), all indexed by the same patient IDs in the same order.

    Returns
    -------
    X_omics : dict[str, pd.DataFrame]
        One patient-by-feature DataFrame per omics view.
    y : pd.Series
        Subtype label per patient, indexed by patient ID.
    """
    omics = pd.read_pickle(data_dir / "omics.pkl")
    y = omics["meta"].astype(str)
    patient_ids = y.index.astype(str)

    X_omics = {name: df for name, df in omics.items() if name != "meta"}

    for name, df in X_omics.items():
        assert df.index.astype(str).equals(patient_ids), f"Patient index differs in {name}"

    return X_omics, y


def make_train_test_split(X_omics, y, test_size, random_state,
                           high_dim_views, n_top_variable_features):
    """Create one shared patient-level train/test split and select variable features.

    High-dimensional views (e.g. transcriptomics, methylation) are reduced to their
    most variable features, computed on the training patients only. Low-dimensional
    views (e.g. proteomics) keep all features. The same selected columns are then
    applied to the test patients.

    Returns
    -------
    X_train_omics, X_test_omics : dict[str, pd.DataFrame]
    y_train, y_test : pd.Series
    train_ids, test_ids : pd.Index
    """
    patient_ids = y.index.astype(str)
    train_ids, test_ids = train_test_split(
        patient_ids, test_size=test_size, random_state=random_state, stratify=y,
    )
    train_ids = pd.Index(train_ids, name="patient_id")
    test_ids = pd.Index(test_ids, name="patient_id")

    X_train_raw = {name: X.loc[train_ids] for name, X in X_omics.items()}
    X_test_raw = {name: X.loc[test_ids] for name, X in X_omics.items()}

    X_train_omics, X_test_omics = {}, {}
    for name, X_train_view in X_train_raw.items():
        if name in high_dim_views:
            n_top = min(n_top_variable_features, X_train_view.shape[1])
            feature_variance = X_train_view.var(axis=0)
            selected_features = feature_variance.sort_values(ascending=False).head(n_top).index
        else:
            selected_features = X_train_view.columns

        X_train_omics[name] = X_train_view.loc[:, selected_features]
        X_test_omics[name] = X_test_raw[name].loc[:, selected_features]

    y_train = y.loc[train_ids]
    y_test = y.loc[test_ids]

    return X_train_omics, X_test_omics, y_train, y_test, train_ids, test_ids


def build_mofa_matrix_input(X_by_omic):
    """Convert aligned omics DataFrames to the nested matrix format expected by mofapy2.

    mofapy2 expects data[view][group] = samples x features. This example uses a
    single group ("TCGA-BRCA_train") and one entry per omics view.
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
    """Fit a compact MOFA model with mofapy2 and return the trained entry point object.

    Uses factor-level ARD so the model can automatically shrink factors that
    explain little signal; the caller later keeps only the "active" factors
    based on the fitted model's variance-explained (R2) table.
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
    """Select 'active' MOFA factors using the fitted model's variance-explained (R2) table.

    Sums each factor's R2 across all omics views and keeps factors whose total R2
    is at least `min_total_r2`. Falls back to the single best factor if none pass.

    Returns
    -------
    active_factor_cols : list[str]
        Factor names (e.g. "Factor1") to keep for downstream analysis.
    factor_r2_summary : pd.DataFrame
        Total R2 per factor, sorted descending.
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


def project_test_patients_to_mofa_factors(model, X_train_by_view, X_test_by_view,
                                           train_factors, view_names):
    """Project held-out patients into the fixed MOFA factor space.

    MOFA is fitted on training patients only. To evaluate on held-out patients we
    keep the learned weights (W) fixed and estimate each test patient's factor
    values via a pseudo-inverse projection, then calibrate the projection to the
    factor scale returned by the trained model. All calibration is fit on
    training patients only.
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
    """Compute one-way ANOVA eta-squared for each factor (fraction of factor variance
    explained by subtype group membership), without fitting a predictive model.
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


def evaluate_predictions(y_true, y_pred, model_name):
    """Return common classification metrics (accuracy, balanced accuracy, macro-F1) in one row."""
    return {
        "model": model_name,
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "macro_f1": f1_score(y_true, y_pred, average="macro"),
    }


def fit_factor_classifier(factors_df, y, train_ids, test_ids, model_name):
    """Fit a logistic regression classifier on training factors and evaluate on
    held-out test factors. Uses only the MOFA factor values (Z), not the raw
    omics matrices, as a diagnostic of whether the learned representation
    preserves subtype information.
    """
    clf = LogisticRegression(
        max_iter=2000,
        class_weight="balanced",
        random_state=RANDOM_STATE,
    )
    clf.fit(factors_df.loc[train_ids], y.loc[train_ids])
    pred = clf.predict(factors_df.loc[test_ids])
    metrics = evaluate_predictions(y.loc[test_ids], pred, model_name)
    return clf, pred, metrics


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_r2_heatmap(r2_all, view_names, active_factor_cols, output_path):
    """Draw and save an annotated heatmap of MOFA variance explained (R2).

    Rows are omics views, columns are the selected "active" factors, and each
    cell is the percentage of that view's variance reconstructed by that
    factor. This is the standard first diagnostic for a fitted MOFA model:
    it shows which views each factor is associated with before inspecting
    patient-level factor values or feature weights.
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
    """Draw and save boxplots of the top subtype-associated MOFA factor values.

    Uses training patients only (matching the eta-squared ranking they were
    selected from). One panel per factor; each panel shows the distribution
    of that factor's values, split out by subtype, to make it easy to see
    which subtypes are shifted relative to one another on that latent axis.
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


def plot_confusion_matrix(y_test, y_pred, output_path):
    """Draw and save a confusion matrix for the held-out subtype predictions.

    Shows how predicted subtypes (from the MOFA-factor logistic regression
    classifier) compare with true subtypes on the test set, making it easy to
    see which subtypes are most often confused with each other.
    """
    fig, ax = plt.subplots(figsize=(5, 5))
    ConfusionMatrixDisplay.from_predictions(y_test, y_pred, xticks_rotation=45, ax=ax)
    ax.set_title("Confusion matrix — MOFA factors")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_ranked_feature_weights(mofa_model_mfx, ranked_factors, view, output_path, n_features=5):
    """Draw and save ranked, signed feature-weight plots for one omics view.

    For each of `ranked_factors`, ranks that view's features by signed MOFA
    weight and plots the curve, labelling the top `n_features` at each tail.
    Features far from zero define the factor most strongly; this makes the
    strongest positive- and negative-contributing features easy to spot for
    each factor.
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


def generate_diagnostic_plots(mofa_model_mfx, factors_df, factor_r2_summary, view_names,
                               active_factor_cols, factor_subtype_assoc, train_ids, y_train,
                               y_test, mofa_pred, output_dir, top_view_for_weights="transcriptomics"):
    """Generate and save the core MOFA diagnostic plots from the notebook.

    This is a single entry point that wraps four complementary diagnostics,
    each saved as its own PNG file under `output_dir`:

    1. R2 heatmap        -- which omics views each active factor reconstructs.
    2. Factor boxplots    -- how the top subtype-associated factors vary by subtype
                              (training patients only).
    3. Confusion matrix   -- held-out subtype prediction errors from the
                              MOFA-factor logistic regression classifier.
    4. Ranked feature weights -- top positive/negative loadings in
                              `top_view_for_weights` for the top subtype-associated
                              factors.

    Together these mirror the notebook's "R2 -> factor values (Z) -> weights (W)"
    interpretation flow, condensed to one representative plot per stage rather
    than every plot variant shown in the notebook.
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
        y_test, mofa_pred,
        output_dir / "part2_mofa_confusion_matrix.png",
    )

    plot_ranked_feature_weights(
        mofa_model_mfx, top_factors[:3], top_view_for_weights,
        output_dir / "part2_mofa_ranked_weights.png",
    )

    print(f"Saved diagnostic plots to: {output_dir}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

    # 1. Load data
    X_omics, y = load_omics_data(DATA_DIR)
    patient_ids = y.index.astype(str)

    # 2. Shared train/test split + variable feature selection for high-dim views
    X_train_omics, X_test_omics, y_train, y_test, train_ids, test_ids = make_train_test_split(
        X_omics, y, TEST_SIZE, RANDOM_STATE, HIGH_DIMENSIONAL_VIEWS, N_TOP_VARIABLE_HIGH_DIM_FEATURES,
    )
    print(f"Training patients: {len(train_ids)} | Test patients: {len(test_ids)}")

    # 3. Build MOFA input and fit (or reload) the model
    mofa_data, view_names, feature_names, sample_names, group_names = build_mofa_matrix_input(X_train_omics)

    if LOAD_MOFA_HDF5 and MOFA_HDF5_FILE.exists():
        mofa_model_mfx = mfx.mofa_model(str(MOFA_HDF5_FILE))
        print(f"Loaded MOFA model from: {MOFA_HDF5_FILE}")
    else:
        fit_mofa(
            mofa_data, view_names, feature_names, sample_names, group_names,
            MAX_FACTORS, MOFA_ITERATIONS, RANDOM_STATE, MOFA_HDF5_FILE,
        )
        mofa_model_mfx = mfx.mofa_model(str(MOFA_HDF5_FILE))
        print(f"Fitted and saved MOFA model to: {MOFA_HDF5_FILE}")

    # 4. Select active factors from the R2 table
    active_factor_cols, factor_r2_summary = select_active_factors(
        mofa_model_mfx, MIN_TOTAL_R2, MAX_FACTORS,
    )
    print(f"Selected {len(active_factor_cols)} active factors out of {MAX_FACTORS} initialized factors.")

    # 5. Read training factors and project test patients into the same factor space
    train_factors_mfx = mofa_model_mfx.get_factors(df=True)
    train_factors_mfx.index = train_factors_mfx.index.astype(str)

    test_factors_mfx = project_test_patients_to_mofa_factors(
        mofa_model_mfx, X_train_omics, X_test_omics, train_factors_mfx, view_names,
    )

    factors_df = pd.concat([train_factors_mfx, test_factors_mfx], axis=0)
    factors_df.columns = factors_df.columns.astype(str)
    factors_df = factors_df.reindex(patient_ids)

    feature_cols = [f for f in active_factor_cols if f in factors_df.columns] or factors_df.columns.tolist()
    print(f"Using {len(feature_cols)} active MOFA factors downstream: {feature_cols}")

    # 6. Quantify factor <-> subtype association (training patients only)
    factor_subtype_assoc = eta_squared_by_factor(
        factors_df.loc[train_ids.astype(str), feature_cols], y_train,
    )
    print("\nTop subtype-associated factors:")
    print(factor_subtype_assoc.head())

    # 7. Predict subtype from MOFA factors
    _, mofa_pred, mofa_metrics = fit_factor_classifier(
        factors_df[feature_cols], y, train_ids.astype(str), test_ids.astype(str),
        "MOFA factors + logistic regression",
    )
    results_df = pd.DataFrame([mofa_metrics])
    print("\nTest-set performance:")
    print(results_df)

    # 8. Save outputs
    predictions_df = pd.DataFrame(
        {
            "patient_id": test_ids.astype(str),
            "true_subtype": y_test.values,
            "mofa_factor_prediction": mofa_pred,
        }
    ).set_index("patient_id")

    factors_to_save = factors_df.copy()
    factors_to_save["split"] = "train"
    factors_to_save.loc[test_ids.astype(str), "split"] = "test"

    results_df.to_csv(OUTPUT_DIR / "part2_mofa_metrics.csv", index=False)
    predictions_df.to_csv(OUTPUT_DIR / "part2_mofa_predictions.csv")
    factors_to_save.to_csv(OUTPUT_DIR / "part2_mofa_factors.csv")
    factor_subtype_assoc.to_csv(OUTPUT_DIR / "part2_mofa_factor_subtype_associations.csv", index=False)

    print("\nSaved outputs to:", OUTPUT_DIR)

    generate_diagnostic_plots(mofa_model_mfx, factors_df, factor_r2_summary, view_names,
                               active_factor_cols, factor_subtype_assoc, train_ids, y_train,
                               y_test, mofa_pred, OUTPUT_DIR, top_view_for_weights="transcriptomics")


if __name__ == "__main__":
    main()