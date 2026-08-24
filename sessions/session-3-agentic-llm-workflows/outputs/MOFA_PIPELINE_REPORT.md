# MOFA Breast Cancer Subtype Prediction Pipeline — Full Results Report

## Executive Summary

The MOFA (Multi-Omics Factor Analysis) pipeline was successfully executed on a cached, pre-trained model using 603 TCGA breast cancer patients with known PAM50 molecular subtypes. The pipeline integrated transcriptomics (2,000 features), proteomics (464 features), and methylation (2,000 features) data, projecting them into a compact 10-factor latent space for subtype classification on held-out test patients.

---

## Key Results

### 1. Active Factors Identified

**All 10 initialized factors are active** (total variance explained ≥ 2.5%):

| Factor | Total R² (%) | Interpretation |
|--------|-------------|-----------------|
| Factor1 | 38.1 | Dominant general variance driver (esp. methylation: 31.2%) |
| Factor2 | 23.1 | **Primary subtype discriminator** (eta² = 0.74) |
| Factor3 | 16.9 | Secondary general variance (balanced across views) |
| Factor4 | 11.4 | Proteomics-associated (16.8% of proteomics R²) |
| Factor5-10 | 3.0-8.6 | Specialized signals in individual views |

**Key observation**: Factor2 is remarkably enriched for subtype signal despite explaining only 23% of total variance, indicating that subtype-associated variation is concentrated in a single latent dimension.

---

### 2. Most PAM50-Associated Factor

**Factor2 is the dominant subtype discriminator:**
- **Eta-squared (η²) = 0.7426** → 74.3% of Factor2's variance is explained by subtype membership
- This is ~2.1× stronger than Factor1 (η² = 0.3498)
- Biological interpretation: Factor2 likely captures hormone receptor/HER2 signaling axes that define PAM50 subtypes

**Ranking of all factors by subtype association:**
```
Factor2 (0.7426) >> Factor1 (0.3498) > Factor7 (0.2642) > Factor8 (0.1732) > ...
```

---

### 3. Held-Out Test Set Classification Performance

**Dataset**: 151 held-out test patients (25% of 603 total, stratified by subtype)

| Metric | Value |
|--------|-------|
| **Overall Accuracy** | **82.8%** |
| **Balanced Accuracy** | **74.1%** |
| **Macro-F1** | **74.7%** |

**Per-subtype performance:**

| Subtype | Correct | Total | Accuracy | Notes |
|---------|---------|-------|----------|-------|
| **Basal** | 22/24 | 91.7% | Clinically distinct; strong signal |
| **HER2** | 9/10 | 90.0% | Well-separated from other subtypes |
| **LumA** | 72/81 | 88.9% | Largest test set; good separation |
| **LumB** | 20/30 | 66.7% | Heterogeneous group; frequent confusion with LumA |
| **Normal** | 2/6 | 33.3% | Very small test set (N=6); background noise in cancer cohort |

**Confusion pattern**: 
- LumB samples are frequently misclassified as LumA (8 instances), reflecting their biological overlap
- Normal samples mostly confused with LumA/LumB (likely due to luminal expression signatures in low-purity samples)
- Basal, HER2, and LumA are cleanly separated

---

## Interpretation

The MOFA model successfully captures breast cancer subtype heterogeneity across transcriptomics, proteomics, and methylation using a minimal 10-factor latent space, with Factor2 emerging as the dominant molecular axis differentiating PAM50 subtypes (η²=0.74). The model achieves 82.8% accuracy on held-out test patients, demonstrating that the learned latent factors preserve strong subtype information despite aggressive dimensionality reduction from >4,464 to 10 features. Performance is particularly robust for luminal A (88.9%), basal (91.7%), and HER2+ (90.0%) subtypes, which are clinically and molecularly well-defined; lower performance on luminal B and normal-like subtypes reflects the known biological heterogeneity within these groups and the inherent difficulty in distinguishing normal samples in a cancer-enriched cohort. The strong association of Factor2 with subtype classification, combined with high balanced accuracy (74.1%) despite class imbalance, suggests that MOFA has identified genuinely interpretable drivers of breast cancer molecular diversity that could serve as targets for subtype-specific stratification and therapy.

---

## Output Files Generated

| File | Contents |
|------|----------|
| `part2_mofa_metrics.csv` | Summary accuracy, balanced accuracy, F1 |
| `part2_mofa_predictions.csv` | Per-patient predictions on test set |
| `part2_mofa_factors.csv` | All factor values for train/test split |
| `part2_mofa_factor_subtype_associations.csv` | Eta-squared values ranking factors |
| `part2_mofa_r2_heatmap.png` | Variance explained by each factor in each view |
| `part2_mofa_confusion_matrix.png` | Test set prediction errors by subtype |
| `part2_mofa_factor_boxplots.png` | Top 4 subtype-associated factors' distributions |
| `part2_mofa_ranked_weights.png` | Top gene/feature loadings for Factors 1,2,7 |

---

## Diagnostic Insights from Plots

1. **R² Heatmap**: Factor1 dominates methylation (31.2%), Factor2 balanced across all views, Factor3 strong in transcriptomics (10.3%)
2. **Confusion Matrix**: Diagonal dominance except LumB↔LumA confusion, validating model's subtype discrimination
3. **Factor Boxplots**: Factor2 shows clear separation by subtype (especially Basal vs. Luminal), supporting high η²
4. **Ranked Weights**: Factor2's top genes include known PAM50 markers; Factor1 methylation-centric; Factor7 transcription-specific

---

## Recommendations for Further Investigation

1. **Investigate Factor2 loadings**: Identify top transcriptomic, proteomic, and methylation markers driving subtype discrimination
2. **Subtype refinement**: Consider re-stratifying LumB/LumA overlap using Factor2 threshold rather than clinical criteria
3. **Normal-like samples**: Filter low-purity or non-malignant samples before retraining; or use Factor2-based confidence scores
4. **Validation**: Test model on independent cohorts (e.g., METABRIC, SCAN-B) to assess generalization

