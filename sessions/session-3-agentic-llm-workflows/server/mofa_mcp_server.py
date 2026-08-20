"""MCP server for the ECCB 2026 MOFA multi-omics practical.

Run with:

    python server/mofa_mcp_server.py

Exposes the same analysis tools as notebook 01, but over MCP (stdio transport).
The heavy work — loading the aligned omics and the pre-fitted MOFA model — is
done ONCE at server startup; the tools are cheap, read-only analyses of that
already-fitted model. Fitting MOFA is deliberately NOT exposed as a tool.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import mofax as mfx

from mcp.server.fastmcp import FastMCP
from mcp.types import ToolAnnotations

# Make `src` importable regardless of the client's working directory.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.mofa_tools import (  # noqa: E402
    RANDOM_STATE, TEST_SIZE, N_TOP_VARIABLE_HIGH_DIM_FEATURES,
    HIGH_DIMENSIONAL_VIEWS, MAX_FACTORS, MIN_TOTAL_R2,
    load_omics_data, make_train_test_split, select_active_factors,
    project_test_patients_to_mofa_factors, eta_squared_by_factor,
    fit_factor_classifier,
)

# ---------------------------------------------------------------------------
# One-time setup: load data + the cached MOFA model (never re-fit).
# ---------------------------------------------------------------------------
DATA_DIR = PROJECT_ROOT / "data"
CACHE = PROJECT_ROOT / "outputs" / (
    f"trained_mofaplus_train_var{N_TOP_VARIABLE_HIGH_DIM_FEATURES}"
    f"_max{MAX_FACTORS}_ard_model.hdf5")

X_OMICS, Y = load_omics_data(DATA_DIR)
X_TRAIN, X_TEST, Y_TRAIN, Y_TEST, TRAIN_IDS, TEST_IDS = make_train_test_split(
    X_OMICS, Y, TEST_SIZE, RANDOM_STATE, HIGH_DIMENSIONAL_VIEWS,
    N_TOP_VARIABLE_HIGH_DIM_FEATURES)
VIEW_NAMES = list(X_TRAIN.keys())
TRAIN_STR, TEST_STR = TRAIN_IDS.astype(str), TEST_IDS.astype(str)

MODEL = mfx.mofa_model(str(CACHE))
_train_factors = MODEL.get_factors(df=True)
_train_factors.index = _train_factors.index.astype(str)
_test_factors = project_test_patients_to_mofa_factors(
    MODEL, X_TRAIN, X_TEST, _train_factors, VIEW_NAMES)
FACTORS = pd.concat([_train_factors, _test_factors], axis=0)
FACTORS.columns = FACTORS.columns.astype(str)
FACTORS = FACTORS.reindex(Y.index.astype(str))
ACTIVE_COLS, R2_SUMMARY = select_active_factors(MODEL, MIN_TOTAL_R2, MAX_FACTORS)
R2_ALL = MODEL.get_r2().rename(columns={"Factor": "factor", "View": "view",
                                        "Group": "group", "R2": "r2"})

mcp = FastMCP("eccb2026-mofa")

# Every tool only READS the fitted model — advertise that to clients.
READ_ONLY = ToolAnnotations(readOnlyHint=True, destructiveHint=False, idempotentHint=True)


# ---------------------------------------------------------------------------
# RESOURCE: a compact summary of the fitted model, to read into context.
# ---------------------------------------------------------------------------
@mcp.resource("mofa://summary")
def mofa_summary() -> dict[str, Any]:
    """Compact summary of the fitted MOFA model and cohort."""
    return {
        "n_patients": int(len(Y)),
        "views": VIEW_NAMES,
        "n_factors": int(MAX_FACTORS),
        "active_factors": ACTIVE_COLS,
        "subtypes": Y.value_counts().to_dict(),
    }


# ---------------------------------------------------------------------------
# PROMPT: a reusable template for interpreting one factor.
# ---------------------------------------------------------------------------
@mcp.prompt()
def interpret_factor(factor: str) -> str:
    """Prompt template asking the model to interpret one MOFA factor end-to-end."""
    return (
        f"Interpret MOFA {factor}. Report which omics views it explains "
        f"(factor_view_r2), how strongly it associates with PAM50 subtype "
        f"(factor_subtype_association), and its strongest transcriptomic drivers "
        f"(top_features_for_factor). Ground every claim in tool output."
    )


# ---------------------------------------------------------------------------
# TOOLS: cheap, read-only analyses of the already-fitted model.
# ---------------------------------------------------------------------------
@mcp.tool(annotations=READ_ONLY)
def data_summary() -> dict:
    """Patients and features per omics view, plus PAM50 subtype class counts."""
    return {"views": {v: {"patients": int(X_OMICS[v].shape[0]),
                          "features": int(X_OMICS[v].shape[1])} for v in VIEW_NAMES},
            "subtype_counts": Y.value_counts().to_dict()}


@mcp.tool(annotations=READ_ONLY)
def split_summary() -> dict:
    """Train/test patient counts and features kept per view after selection."""
    return {"n_train": int(len(TRAIN_IDS)), "n_test": int(len(TEST_IDS)),
            "features_per_view": {v: int(X_TRAIN[v].shape[1]) for v in VIEW_NAMES}}


@mcp.tool(annotations=READ_ONLY)
def active_factors() -> dict:
    """Which MOFA factors are active and each factor's total R2 across views."""
    return {"active_factors": ACTIVE_COLS,
            "total_r2": R2_SUMMARY.set_index("factor")["total_r2"].round(3).to_dict()}


@mcp.tool(annotations=READ_ONLY)
def factor_view_r2(factor: str) -> dict:
    """Variance explained (R2) by one factor in each omics view, and its top view."""
    sub = R2_ALL[R2_ALL["factor"] == factor][["view", "r2"]]
    d = {k: round(float(v), 3) for k, v in zip(sub["view"], sub["r2"])}
    return {"factor": factor, "r2_by_view": d,
            "top_view": max(d, key=d.get) if d else None}


@mcp.tool(annotations=READ_ONLY)
def factor_subtype_association() -> list:
    """Rank active factors by eta-squared association with PAM50 subtype (train)."""
    assoc = eta_squared_by_factor(FACTORS.loc[TRAIN_STR, ACTIVE_COLS], Y_TRAIN)
    return assoc.assign(eta_squared=assoc["eta_squared"].round(3)).to_dict("records")


@mcp.tool(annotations=READ_ONLY)
def top_features_for_factor(factor: str, view: str = "transcriptomics", n: int = 5) -> dict:
    """Top positive/negative weighted features of a view for a factor (its drivers)."""
    w = MODEL.get_weights(views=view, df=True)
    w.columns = w.columns.astype(str)
    s = w[factor].sort_values()
    return {"factor": factor, "view": view,
            "top_negative": {k: round(float(v), 3) for k, v in s.head(n).items()},
            "top_positive": {k: round(float(v), 3) for k, v in s.tail(n).items()}}


@mcp.tool(annotations=READ_ONLY)
def classify_subtype_from_factors() -> dict:
    """Predict subtype from MOFA factors; held-out metrics + most-confused pair."""
    _, pred, metrics = fit_factor_classifier(
        FACTORS[ACTIVE_COLS], Y, TRAIN_STR, TEST_STR, "MOFA factors + LR")
    labels = sorted(set(Y_TEST.values) | set(pred))
    cm = pd.crosstab(pd.Series(Y_TEST.values, name="true"),
                     pd.Series(pred, name="pred")).reindex(
        index=labels, columns=labels, fill_value=0)
    off = cm.to_numpy(copy=True)
    np.fill_diagonal(off, 0)
    r, c = np.unravel_index(off.argmax(), off.shape)
    return {"metrics": {k: round(v, 3) for k, v in metrics.items() if k != "model"},
            "most_confused": {"true": labels[r], "predicted": labels[c], "count": int(off.max())}}


@mcp.tool(annotations=READ_ONLY)
def train_vs_test_subtype_association() -> list:
    """Compare factor<->subtype eta-squared on train vs projected test patients."""
    tr = eta_squared_by_factor(FACTORS.loc[TRAIN_STR, ACTIVE_COLS], Y_TRAIN)
    te = eta_squared_by_factor(FACTORS.loc[TEST_STR, ACTIVE_COLS], Y_TEST)
    m = tr.merge(te, on="factor", suffixes=("_train", "_test"))
    m["eta_squared_train"] = m["eta_squared_train"].round(3)
    m["eta_squared_test"] = m["eta_squared_test"].round(3)
    return m.to_dict("records")


if __name__ == "__main__":
    mcp.run()
