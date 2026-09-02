#!/usr/bin/env python
"""
Run the full breast cancer subtype prediction pipeline on the fitted MOFA model.

This script:
1. Loads the cached MOFA model
2. Reports which factors are active
3. Identifies which factor is most associated with PAM50 subtype
4. Reports held-out classification performance
"""

from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

import numpy as np
import pandas as pd
import mofax as mfx

from mofa_tools import (
    load_omics_data,
    make_train_test_split,
    select_active_factors,
    project_test_patients_to_mofa_factors,
    eta_squared_by_factor,
    fit_factor_classifier,
    generate_diagnostic_plots,
)

# Configuration
RANDOM_STATE = 42
SESSION_DIR = Path(__file__).resolve().parent
DATA_DIR = Path("/data") if Path("/data").exists() else SESSION_DIR / "data"
OUTPUT_DIR = SESSION_DIR / "outputs"

TEST_SIZE = 0.25
N_TOP_VARIABLE_HIGH_DIM_FEATURES = 2000
HIGH_DIMENSIONAL_VIEWS = ["transcriptomics", "methylation", "proteomics"]
MAX_FACTORS = 10
MIN_TOTAL_R2 = 2.5
MOFA_HDF5_FILE = OUTPUT_DIR / (
    f"trained_mofaplus_train_var{N_TOP_VARIABLE_HIGH_DIM_FEATURES}"
    f"_max{MAX_FACTORS}_ard_model.hdf5"
)


def main():
    print("=" * 80)
    print("MOFA Breast Cancer Subtype Prediction Pipeline")
    print("=" * 80)
    
    # 1. Load data
    print("\n[1/5] Loading omics data...")
    X_omics, y = load_omics_data(DATA_DIR)
    patient_ids = y.index.astype(str)
    print(f"    Loaded {len(patient_ids)} patients with PAM50 subtypes: {y.unique()}")
    
    # 2. Create shared train/test split
    print(f"\n[2/5] Creating train/test split (test_size={TEST_SIZE})...")
    X_train_omics, X_test_omics, y_train, y_test, train_ids, test_ids = make_train_test_split(
        X_omics, y, TEST_SIZE, RANDOM_STATE, HIGH_DIMENSIONAL_VIEWS, N_TOP_VARIABLE_HIGH_DIM_FEATURES,
    )
    print(f"    Training patients: {len(train_ids)} | Test patients: {len(test_ids)}")
    print(f"    Features selected per view (top {N_TOP_VARIABLE_HIGH_DIM_FEATURES} variable):")
    for name in HIGH_DIMENSIONAL_VIEWS + [v for v in X_train_omics.keys() if v not in HIGH_DIMENSIONAL_VIEWS]:
        print(f"        {name}: {X_train_omics[name].shape[1]} features")
    
    # 3. Load cached MOFA model
    print(f"\n[3/5] Loading cached MOFA model from: {MOFA_HDF5_FILE}")
    if not MOFA_HDF5_FILE.exists():
        print(f"    ERROR: Model file not found at {MOFA_HDF5_FILE}")
        sys.exit(1)
    
    mofa_model_mfx = mfx.mofa_model(str(MOFA_HDF5_FILE))
    # Extract view names from R2 table
    r2_temp = mofa_model_mfx.get_r2().rename(
        columns={"Factor": "factor", "View": "view", "Group": "group_mofax", "R2": "r2"}
    )
    view_names = r2_temp["view"].unique().tolist()
    print(f"    Model loaded with {len(view_names)} views: {view_names}")
    
    # 4. Select active factors
    print(f"\n[4/5] Identifying active factors (min_total_r2={MIN_TOTAL_R2})...")
    active_factor_cols, factor_r2_summary = select_active_factors(
        mofa_model_mfx, MIN_TOTAL_R2, MAX_FACTORS,
    )
    print(f"    Found {len(active_factor_cols)} active factors:")
    print(factor_r2_summary.to_string(index=False))
    
    # 5. Extract training factors and project test patients
    print(f"\n[5/5] Extracting factor values and projecting held-out test patients...")
    train_factors_mfx = mofa_model_mfx.get_factors(df=True)
    train_factors_mfx.index = train_factors_mfx.index.astype(str)
    print(f"    Training factors shape: {train_factors_mfx.shape}")
    
    test_factors_mfx = project_test_patients_to_mofa_factors(
        mofa_model_mfx, X_train_omics, X_test_omics, train_factors_mfx, view_names,
    )
    print(f"    Test factors shape: {test_factors_mfx.shape}")
    
    # Combine train and test factors
    factors_df = pd.concat([train_factors_mfx, test_factors_mfx], axis=0)
    factors_df.columns = factors_df.columns.astype(str)
    factors_df = factors_df.reindex(patient_ids)
    
    feature_cols = [f for f in active_factor_cols if f in factors_df.columns] or factors_df.columns.tolist()
    print(f"    Using {len(feature_cols)} active factors: {feature_cols}")
    
    # =========================================================================
    # MAIN RESULTS REPORTING
    # =========================================================================
    
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)
    
    # Active factors
    print("\n1. ACTIVE FACTORS:")
    print(f"   Total active factors: {len(active_factor_cols)}")
    print(f"   Factors: {', '.join(active_factor_cols)}")
    
    # Factor-subtype associations
    print("\n2. FACTOR-SUBTYPE ASSOCIATIONS (Eta-Squared):")
    print("   Computed on training patients only.")
    factor_subtype_assoc = eta_squared_by_factor(
        factors_df.loc[train_ids.astype(str), feature_cols], y_train,
    )
    print("\n   All active factors ranked by association strength:")
    print(factor_subtype_assoc.to_string(index=False))
    
    most_associated = factor_subtype_assoc.iloc[0]
    print(f"\n   *** Most associated factor: {most_associated['factor']} (eta² = {most_associated['eta_squared']:.4f})")
    
    # Classification performance
    print("\n3. HELD-OUT TEST SET CLASSIFICATION PERFORMANCE:")
    print("   Logistic Regression classifier trained on MOFA factors")
    _, mofa_pred, mofa_metrics = fit_factor_classifier(
        factors_df[feature_cols], y, train_ids.astype(str), test_ids.astype(str),
        "MOFA factors + logistic regression",
    )
    
    print(f"   Accuracy:           {mofa_metrics['accuracy']:.4f}")
    print(f"   Balanced Accuracy:  {mofa_metrics['balanced_accuracy']:.4f}")
    print(f"   Macro-averaged F1:  {mofa_metrics['macro_f1']:.4f}")
    
    # Generate and save diagnostic plots
    print("\n4. GENERATING DIAGNOSTIC PLOTS...")
    generate_diagnostic_plots(
        mofa_model_mfx, factors_df, factor_r2_summary, view_names,
        active_factor_cols, factor_subtype_assoc, train_ids, y_train,
        y_test, mofa_pred, OUTPUT_DIR, top_view_for_weights="transcriptomics"
    )
    
    # Save results tables
    print("\n5. SAVING RESULTS TABLES TO OUTPUTS...")
    OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    
    results_df = pd.DataFrame([mofa_metrics])
    results_df.to_csv(OUTPUT_DIR / "part2_mofa_metrics.csv", index=False)
    print(f"    Saved metrics to: part2_mofa_metrics.csv")
    
    predictions_df = pd.DataFrame({
        "patient_id": test_ids.astype(str),
        "true_subtype": y_test.values,
        "mofa_factor_prediction": mofa_pred,
    }).set_index("patient_id")
    predictions_df.to_csv(OUTPUT_DIR / "part2_mofa_predictions.csv")
    print(f"    Saved predictions to: part2_mofa_predictions.csv")
    
    factors_to_save = factors_df.copy()
    factors_to_save["split"] = "train"
    factors_to_save.loc[test_ids.astype(str), "split"] = "test"
    factors_to_save.to_csv(OUTPUT_DIR / "part2_mofa_factors.csv")
    print(f"    Saved factor values to: part2_mofa_factors.csv")
    
    factor_subtype_assoc.to_csv(OUTPUT_DIR / "part2_mofa_factor_subtype_associations.csv", index=False)
    print(f"    Saved factor-subtype associations to: part2_mofa_factor_subtype_associations.csv")
    
    print("\n" + "=" * 80)
    print("PIPELINE COMPLETE")
    print("=" * 80)
    
    return {
        'active_factors': active_factor_cols,
        'most_associated_factor': most_associated['factor'],
        'eta_squared': most_associated['eta_squared'],
        'accuracy': mofa_metrics['accuracy'],
        'balanced_accuracy': mofa_metrics['balanced_accuracy'],
        'macro_f1': mofa_metrics['macro_f1'],
    }


if __name__ == "__main__":
    results = main()
