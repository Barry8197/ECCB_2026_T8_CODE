# Glossary 

# Session 2 : Creating Multi-Omics Profiles

## Global concepts across all parts

| Term | Definition |
| --- | --- |
| **Omic / omics layer** | One molecular data type profiled on the same patients — transcriptomics (RNA-seq), proteomics, or methylation. |
| **Subtype** | The supervised prediction target — the breast cancer subtype label (e.g., Basal, LumA, LumB, Her2, Normal) for each patient in TCGA-BRCA. |
| **Patient ID / index** | The identifier used to align rows (patients) consistently across all omics matrices and the label vector. |
| **Feature** | A single measured variable within an omic (e.g., one gene's expression, one protein level, one methylation probe). |
| --- | --- |
| **Train/test split** | A single, patient-level split shared across all omics views, ensuring a patient's data isn't split across train and test in different views. |
| **Out-of-fold (OOF) probabilities** | Predicted probabilities for training patients generated via cross-validation, so that a meta-learner (stacking) is trained on predictions the base model didn't get to "cheat" on — avoiding information leakage. |
| **Cross-validation (`cross_val_predict`)** | The technique used to generate out-of-fold predictions on the training set for building the stacking meta-learner. |
| **Data leakage** | When information from the test set inadvertently influences model training (e.g. fitting a scaler on all data instead of just training data), leading to overly optimistic performance estimates. |
| **Standardization / scaling** | Rescaling features (via `StandardScaler`) to comparable ranges before fitting; done inside the pipeline so scaling parameters are learned only from training data, not test data. |
| **Label encoding** | Converting categorical class labels (e.g. cancer subtypes) into integer codes that a model can work with. |
| **Held-out test set** | The patients not used in training, used only to evaluate final model performance. |
| --- | --- |
| **Balanced accuracy** | The primary comparison metric — the average of per-class recall, which avoids favoring predictions on larger subtype classes. |
| **Macro-F1** | The unweighted average F1 score across all subtype classes. |
| **Classification report** | Per-class precision, recall, and F1-score summary for a model's predictions. |
| **Confusion matrix** | A table cross-tabulating true vs. predicted subtypes|
| **PCA (Principal Component Analysis)** | Used only as a quick, exploratory 2D visualization of each omic (not fed into the classifiers) to see whether different omics separate patients/subtypes differently. |

## Part 1: Linear Methods for Multi-Omic Integration

### Linear integration strategies

| Term | Meaning in this notebook |
| --- | --- |
| **Early integration (concatenation)** | Combining all omics (transcriptomics, proteomics, methylation) into one wide feature matrix before training a single classifier. Simple and often strong, but increases dimensionality and doesn't distinguish shared, redundant, or omic-specific signal. |
| **Late integration (prediction averaging)** | Training one classifier per omic, then combining their predicted class probabilities (by default, an equal-weighted average) into a final prediction. Keeps each omic's model separate and interpretable. |
| **Stacking (learned combiner / meta-learner)** | A more adaptive form of late integration: instead of a fixed averaging rule, a second ("meta") model is trained on the per-omic predicted probabilities to learn how much each omic should contribute, potentially differently per class. |
| **Single-omic baseline** | A classifier trained on just one omics layer, used as a reference point to judge whether any integration strategy actually improves over using the best individual omic alone. |
| **Fusion rule** | The rule used to combine modality-specific outputs — e.g., fixed equal-weight averaging (late integration) versus a learned rule (stacking). |

## Part 2: Correlation-Based / Factor-Based Multi-Omic Integration with MOFA

### Core MOFA concepts

| Term | Definition |
| --- | --- |
| **MOFA (Multi-Omics Factor Analysis)** | An unsupervised, intermediate-integration method that keeps omics views separate but learns one shared patient-level latent representation across them. |
| **Intermediate integration (latent-factor integration)** | An integration strategy that is neither early integration (concatenating all omics into one table) nor late integration (training separate models per view and combining results at the end); MOFA instead learns shared latent factors directly from multiple views. |
| --- | --- |
| **Sample** | One patient. Each row in the omics matrices is a patient. |
| **View** | One omics layer measured on the same patients (e.g., transcriptomics, proteomics, methylation). |
| **Group** | A cohort or batch of samples; here a single group, `TCGA-BRCA`. |
| **Feature** | One measured variable within a view — a gene, protein, or methylation probe. |
| **Factor** | A hidden ("latent") axis of variation learned by MOFA. May reflect a biological process, technical effect, subtype gradient, or other source of covariance across views. Factors are not automatically biological pathways — they must be interpreted after fitting. |
| --- | --- |
| **Factor value** | The coordinate of one patient on one factor — how strongly that factor is active in that patient. Collected in the factor matrix **Z**. |
| **Weight / loading** | The contribution of one feature to one factor — which molecular features define the factor. Collected in the weight matrix **W** (one per view). |
| **Variance explained (R2)** | A summary (as a %) of how much of a view's variance is reconstructed by a given factor. Read from the fitted model via `mofax`. |
| **Factor activity** | A factor's value in a given patient (or the average factor value within a group, e.g., a subtype), used to describe how "on" or "off" that latent axis is. |
| **Variable feature selection** | Reducing high-dimensional views (transcriptomics, methylation) to their most variable features, computed on training patients only |


### MOFA definitions, structure, and math

| Term | Definition|
| --- | --- |
| **X_view / Y(v)** | One molecular table (e.g., transcriptomics) — the observed omics data for view *v*. |
| **Z (factor matrix)** | Patients × factors matrix of factor values, shared across all views. |
| **W(v) (weight/loading matrix)** | Features × factors matrix for view *v*, giving each feature's contribution to each factor. |
| **Reconstruction** | The approximation of each input matrix as `X_view ≈ Z × W(v) + noise`; conceptually similar to PCA but extended to multiple omics views simultaneously. |
| **Sign ambiguity** | The direction (positive/negative) of a factor and its weights is arbitrary — multiplying both by −1 describes the same pattern. Only relative differences between patients/features matter. |
| **ARD (Automatic Relevance Determination)** | A model mechanism (`ard_factors`, `ard_weights`) that automatically shrinks/down-weights factors or features that explain little signal, letting a compact model discover how many factors are actually needed. |
| **Spike-and-slab prior** | A sparsity-inducing prior (`spikeslab_weights`) that lets individual feature weights be effectively zero, making factors easier to interpret by concentrating on a smaller set of strong features. |
| **Variational inference** | The optimization approach MOFA uses to fit its probabilistic model — searching for factors and weights that best reconstruct the omics data while accounting for noise. |
| **Eta-squared (η²)** | An ANOVA-based statistic estimating the proportion of a factor's total variation explained by subtype-group differences; used to rank factors by subtype association. |
| **Active factors** | The subset of the initialized (upper-bound) factors that are kept for interpretation, selected by summing each factor's R2 across views and keeping those above a chosen cutoff (`MIN_TOTAL_R2`). |
| **Projection** | Estimating factor values for new (e.g., held-out test) patients by keeping the learned MOFA weights fixed and computing where those patients fall on the existing factors, without refitting the model. |

## Part 3: Deep Learning Integration Methods Methods

### Core Deep-Learning Concepts

| Term | Definition |
|------|------------|
| **Neural network** | A model made of layers of connected "neurons" (weighted sums followed by a non-linear activation) that learns to map inputs to outputs by adjusting weights during training. |
| **Multi-layer perceptron (MLP)** | A basic type of neural network consisting of one or more fully-connected (`Linear`) layers stacked with non-linear activations in between. |
| **Multi-modal encoder** | A fusion strategy that trains a separate encoder for each data view first, then combines (concatenates) the resulting embeddings before a shared classifier. |
| **Integrated Gradients** | A gradient-based feature-attribution method that estimates each input feature's contribution to a model's prediction by integrating gradients along a path from a baseline input to the actual input. |
| --- | --- |
| **Encoder** | The part of a neural network that compresses raw input features into a lower-dimensional internal representation (an embedding). |
| **Embedding** | A dense, lower-dimensional numeric vector that represents an input (e.g. a patient's omics profile) in a way that preserves useful structure for downstream tasks. |
| **Linear layer** | A neural network layer that applies a weighted sum plus a bias to its inputs (`y = Wx + b`); the basic building block of an MLP. |
| --- | --- |
| **ReLU (Rectified Linear Unit)** | A non-linear activation function, `f(x) = max(0, x)`, applied after linear layers so the network can learn non-linear relationships. |
| **Dropout** | A regularisation technique that randomly "switches off" a fraction of neurons during training to reduce overfitting. |
| **Batch / mini-batch** | A small subset of training samples processed together in one forward/backward pass, rather than using the whole dataset at once. |
| **Epoch** | One complete pass through the entire training dataset during model training. |
| **DataLoader** | A PyTorch utility that wraps a dataset and serves it in shuffled mini-batches during training/evaluation. |
| **Cross-entropy loss** | A loss function commonly used for multi-class classification; it measures how well the predicted class probabilities match the true labels. |
| **Adam optimiser** | An adaptive gradient-descent algorithm commonly used to update neural network weights during training. |
| **Weight decay (L2 regularisation)** | A penalty added during training that discourages large weight values, helping to reduce overfitting. |
| **Gradient descent** | An optimisation method that iteratively adjusts model parameters in the direction that reduces the loss function. |
