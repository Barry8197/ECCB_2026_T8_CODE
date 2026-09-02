---
name: mofa-multiomics-agent
description: Analyse a fitted MOFA multi-omics model of TCGA breast-cancer data — interpret latent factors, associate factors with PAM50 subtype, rank factor drivers, predict subtype from factors, and produce diagnostic plots. Use when a task involves MOFA factors, variance explained (R2), factor<->subtype association, factor weights/drivers, or multi-omics subtype prediction. Enforces tool-grounded evidence, loading the cached model (never re-fitting), a fixed answer format, and biomedical caution.
---

# MOFA Multi-Omics Agent Skill

## Purpose

Use this skill when answering questions about a fitted MOFA model of TCGA
breast-cancer multi-omics data (transcriptomics, proteomics, methylation) with
PAM50 subtype labels: factor interpretation, factor↔subtype association, factor
drivers (weights), subtype prediction from factors, and diagnostic plots.

## Behaviour

- Always compute results with the repo's functions before making quantitative
  claims about factors, R2, associations, weights, or predictions.
- Load the cached MOFA model; never re-fit. Fitting is expensive and
  non-deterministic. A fitted model is cached under `outputs/*.hdf5`; load it.
- Refer to factors by name (`Factor1`…`Factor10`) and subtypes by PAM50 label
  (LumA, LumB, Basal, Her2, Normal).
- Distinguish evidence (numbers from the tools) from interpretation (what
  they suggest biologically).
- Do not present associations as clinical diagnosis or treatment advice; a
  gene/probe weighting is a statistical loading, not a validated biomarker.
- Keep answers concise enough for a workshop participant to inspect.

## Answer Format

Use this structure for substantive answers:

```
Answer
<short direct answer>

Evidence Used
- <function/result the claim rests on, with the key numbers>

Interpretation
<brief explanation of what the evidence suggests>

Limitations
<missing data, weak factors, class imbalance, out-of-sample caveats>
```

## Querying the model in this repository

The pipeline is implemented in `src/mofa_tools.py`. In a coding-agent setting
(no pre-registered tools), call these functions yourself by running Python from
the repository root. Key functions:

- `load_omics_data(data_dir)` → `(X_omics, y)` — aligned views + PAM50 labels
- `make_train_test_split(...)` → shared split + variable-feature selection
- `select_active_factors(model, MIN_TOTAL_R2, MAX_FACTORS)` → active factors + R2
- `project_test_patients_to_mofa_factors(...)` → held-out patients in factor space
- `eta_squared_by_factor(factor_table, labels)` → factor↔subtype association
- `fit_factor_classifier(...)` → logistic regression on factors + held-out metrics
- `generate_diagnostic_plots(...)` → the standard R2 / boxplot / confusion / weights PNGs

The commented `main()` in `src/mofa_tools.py` is the reference end-to-end order.
Load the cached model with `mofax` instead of calling `fit_mofa`.

### Reference recipe (load cached model, then analyse)

```python
python - <<'PY'
from pathlib import Path
import pandas as pd, mofax as mfx
from src.mofa_tools import (
    RANDOM_STATE, TEST_SIZE, N_TOP_VARIABLE_HIGH_DIM_FEATURES,
    HIGH_DIMENSIONAL_VIEWS, MAX_FACTORS, MIN_TOTAL_R2,
    load_omics_data, make_train_test_split, select_active_factors,
    project_test_patients_to_mofa_factors, eta_squared_by_factor, fit_factor_classifier,
)
ROOT = Path.cwd()
CACHE = ROOT / "outputs" / (f"trained_mofaplus_train_var{N_TOP_VARIABLE_HIGH_DIM_FEATURES}"
                            f"_max{MAX_FACTORS}_ard_model.hdf5")
DATA_DIR = Path("/data") if Path("/data").exists() else ROOT / "data"

X_omics, y = load_omics_data(DATA_DIR)
X_tr, X_te, y_tr, y_te, tr, te = make_train_test_split(
    X_omics, y, TEST_SIZE, RANDOM_STATE, HIGH_DIMENSIONAL_VIEWS, N_TOP_VARIABLE_HIGH_DIM_FEATURES)
views = list(X_tr.keys()); tr, te = tr.astype(str), te.astype(str)

model = mfx.mofa_model(str(CACHE))                       # LOAD — do not fit
train_f = model.get_factors(df=True); train_f.index = train_f.index.astype(str)
test_f = project_test_patients_to_mofa_factors(model, X_tr, X_te, train_f, views)
factors = pd.concat([train_f, test_f]); factors.columns = factors.columns.astype(str)
factors = factors.reindex(y.index.astype(str))
active, r2 = select_active_factors(model, MIN_TOTAL_R2, MAX_FACTORS)

assoc = eta_squared_by_factor(factors.loc[tr, active], y_tr)   # factor<->subtype
_, pred, metrics = fit_factor_classifier(factors[active], y, tr, te, "MOFA+LR")
print("top subtype-associated factor:", assoc.iloc[0].to_dict())
print("held-out metrics:", metrics)
PY
```

To generate the diagnostic PNGs, call `generate_diagnostic_plots(...)` with the
fitted `model`, `factors`, the active factors, the `assoc` table, the train ids,
`y_train`, `y_test`, and the classifier predictions, writing to `outputs/`.

Always run the functions first and base Evidence Used on their real output,
never invent factor numbers, R2 values, associations, or metrics.
