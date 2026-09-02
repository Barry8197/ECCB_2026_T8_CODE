#!/usr/bin/env python
"""One-time (re)fit of the cached MOFA model against the local data.

Run this whenever omics.pkl changes so the cached hdf5 in outputs/ matches the
data participants load (the notebooks and MCP server assert the cache exists
and never re-fit). Uses the same constants and split as the notebooks, and the
same DATA_DIR resolution (/data if present, else the clone's own data/).

Usage:
    python retrain_mofa_cache.py            # fit and overwrite outputs/<cache>.hdf5
    python retrain_mofa_cache.py --dry-run  # load + split + report only, no fit
"""
from pathlib import Path
import sys

SESSION_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SESSION_DIR))

from src.mofa_tools import (
    RANDOM_STATE, TEST_SIZE, N_TOP_VARIABLE_HIGH_DIM_FEATURES,
    HIGH_DIMENSIONAL_VIEWS, MAX_FACTORS,
    load_omics_data, make_train_test_split, build_mofa_matrix_input, fit_mofa,
)

DATA_DIR = Path("/data") if Path("/data").exists() else SESSION_DIR / "data"
OUTFILE = SESSION_DIR / "outputs" / (
    f"trained_mofaplus_train_var{N_TOP_VARIABLE_HIGH_DIM_FEATURES}"
    f"_max{MAX_FACTORS}_ard_model.hdf5")
ITERATIONS = 500  # matches the original cache's training_opts (converges long before)


def main(dry_run: bool) -> None:
    print(f"data dir : {DATA_DIR}")
    X_omics, y = load_omics_data(DATA_DIR)
    print(f"patients : {len(y)}  subtypes: {sorted(y.unique())}")

    X_train, X_test, y_train, y_test, train_ids, test_ids = make_train_test_split(
        X_omics, y, TEST_SIZE, RANDOM_STATE, HIGH_DIMENSIONAL_VIEWS,
        N_TOP_VARIABLE_HIGH_DIM_FEATURES)
    print(f"split    : {len(train_ids)} train / {len(test_ids)} test")
    print(f"views    : {[f'{k}:{v.shape[1]}' for k, v in X_train.items()]}")

    if dry_run:
        print("dry run — not fitting.")
        return

    data, view_names, feature_names, sample_names, group_names = \
        build_mofa_matrix_input(X_train)
    if OUTFILE.exists():
        OUTFILE.unlink()  # fit_mofa only saves when the outfile is absent
    print(f"fitting MOFA (max {MAX_FACTORS} factors, {ITERATIONS} max iters) ...")
    fit_mofa(data, view_names, feature_names, sample_names, group_names,
             MAX_FACTORS, ITERATIONS, RANDOM_STATE, OUTFILE)

    import mofax as mfx
    m = mfx.mofa_model(str(OUTFILE))
    n_model = m.get_factors(df=True).shape[0]
    m.close()
    assert n_model == len(train_ids), (
        f"saved model has {n_model} samples but the split has {len(train_ids)} - save failed?")
    print(f"saved    : {OUTFILE}  ({n_model} training samples, verified)")


if __name__ == "__main__":
    main(dry_run="--dry-run" in sys.argv)
