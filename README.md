# Identifying Hypoxia-Responsive Genes from Breast Cancer RNA-seq Data

## Overview

As a part of my senior year internship, I made this project to investigate transcriptional responses to hypoxia in breast cancer cell lines using RNA-seq data from GEO accession GSE111653. This project has helped me practice various methods of analysis used in simple bioinformatics processes.

The goal was to:

1. Identify genes differentially expressed under hypoxic conditions.
2. Visualize transcriptomic variation using PCA.
3. Build a machine-learning model capable of distinguishing hypoxic and normoxic samples.
4. Evaluate whether a compact gene signature can classify hypoxic and normoxic samples within this dataset.
---

## Dataset

Source:

GEO: GSE111653

"A Consensus Hypoxia Signature in Breast Cancer"

The dataset contains RNA-seq expression profiles from breast cancer cell lines cultured under:

- Normoxia (20% O2)
- Hypoxia (1% O2)

for 24 hours.

---

## Initial Methods

### Preprocessing

- Log2 transformation (under pca_analysis)
- Low-expression filtering

### Exploratory Analysis

- Principal Component Analysis (PCA)

### Differential Expression

Mean expression differences between hypoxic and normoxic samples were calculated to identify candidate hypoxia-responsive genes.

### Machine Learning

Random Forest classifiers were trained using:

- All genes (file not included in the repository due to low accuracy level - high noise)
- Top differentially expressed genes (classifier_top_genes.py - this file had a minor data leakage, see updated methodology section below)

Model performance was evaluated using:

- Accuracy
- ROC-AUC
- 5-fold cross-validation

---

## Key Results

### Differential Expression

Strongly induced hypoxia genes included:

- CA9
- NDRG1
- EGLN3
- BNIP3
- P4HA1
- ANGPTL4
- ADM
- LOX

### Classification Performance

| Feature Set | Accuracy |
|-------------|----------|  
| All genes   |   38.5%  |
| Top 200     |   76.9%  |
| Top 100     |   84.6%  |
| Top 50      |   84.6%  |
| Top 20      |   92.3%  |

### Cross-Validation

- Mean Accuracy: 82.4%
- Mean ROC-AUC: 0.863

### Feature Importance

Top predictive genes:

1. ADM
2. P4HA1
3. BNIP3
4. CASP14
5. TCAF2

---
## Main Findings

1. Several well-known hypoxia-related genes, including CA9, NDRG1, BNIP3, P4HA1, EGLN3, ANGPTL4, LOX, and ADM, showed much higher expression under hypoxic conditions than under normoxic conditions.

2. PCA analysis suggested that differences between breast cancer cell lines had a larger impact on overall gene expression patterns than oxygen levels alone.

3. A machine-learning model trained on all ~56,000 genes performed poorly (38.5% accuracy). However, limiting the model to the top 20 hypoxia-responsive genes greatly improved performance.

4. Using five-fold cross-validation, I found that the final model achieved an average accuracy of 82.4% and an average ROC-AUC score of 0.863, indicating that a small set of genes can generally distinguish hypoxic and normoxic samples within this dataset.

## Interpretation

A small set of hypoxia-responsive genes substantially outperformed models trained on the full transcriptome.

This can suggest that a compact gene signature captures most of the predictive information associated with exposure to a hypoxic environment.

The identified genes are consistent with previously reported hypoxia-regulated pathways in breast cancer (https://pubmed.ncbi.nlm.nih.gov/24156323/ was one of my sources).

---

## Limitations

- Small sample size (62 samples)
- No external validation cohort
- Cell-line models may not fully represent patient tumors
- Feature selection was performed before cross-validation
- Data leakage: potentially inflated results

---

## Methodology Note: Identifying and Fixing Data Leakage

I noticed that the original classification pipeline (`classifier_top_genes.py` and `cross_validation.py`)
had two methodological issues that can inflate reported performance:

1. **Gene selection used the whole dataset, including the test data**: Before splitting
   the data into train/test groups, I picked the top genes using *all* the samples,
   including the ones that were later used to test the model. This means the test data
   wasn't really "unseen" by the time evaluation happened, since it helped decide which
   genes got picked in the first place. That can make the model look better than it
   actually is at predicting new data.

2. **Matched samples got split apart by accident**: Each cell line has two samples: one
   grown in hypoxia and one in normoxia. The original split didn't account for this, so
   sometimes one sample from a pair ended up in training while its match ended up in
   testing. This could let the model partly "cheat" by recognizing which cell line a
   sample came from (since it already saw the matching sample during training), instead
   of only learning the actual hypoxia vs. normoxia difference.

To address both issues, `nested_cv_classifier.py` implements **grouped, nested
cross-validation**:
- Samples are split using `StratifiedGroupKFold`, grouped by cell line, so paired
  samples are never separated across train/test.
- Differential expression (gene selection) is recomputed *independently within each
  training fold*, using only that fold's training cell lines — never the held-out
  test cell lines.

**All other analyses in this repository (feature importance, ROC curve, top-gene
plots) were performed using the original `classifier_top_genes.py` pipeline and
predate this correction.** They are kept as-is for transparency, but the headline
performance numbers below come from the corrected, leakage-free pipeline.

### Leaky vs. Fixed Methodology (TOP_N = 20, held constant for fair comparison)

| Metric    | Original (leaky, ungrouped) | Fixed (nested, grouped) |
|-----------|------------------------------|--------------------------|
| Accuracy  | 82.4%                        | 86.0%                    |
| ROC-AUC   | 0.863                        | 0.920                    |


Despite removing both sources of leakage and evaluating on entirely unseen cell
lines, performance didn't drop - it held steady / slightly improved. This
suggests the hypoxia-responsive gene signature genuinely generalizes to new,
unseen cell lines, rather than the original numbers being artifacts of leakage.

### Per-Fold Results (TOP_N = 20)

| Fold | Accuracy | ROC-AUC | Test cell lines (n) |
|------|----------|---------|----------------------|
| 1    | 0.714    | 0.796   | 7                    |
| 2    | 0.917    | 0.917   | 6                    |
| 3    | 0.917    | 1.000   | 6                    |
| 4    | 0.833    | 0.889   | 6                    |
| 5    | 0.917    | 1.000   | 6                    |

**Mean Accuracy: 86.0% (± see script output)**
**Mean ROC-AUC: 0.920 (± see script output)**

Note: perfect scores in folds 3 and 5 reflect small per-fold test sets (n=6 cell
lines) rather than indicating an error-free model — a note worth keeping in
mind given the limited overall sample size.

## References

Ye IC, Fertig EJ, DiGiacomo JW, et al.

Molecular Portrait of Hypoxia in Breast Cancer: A Prognostic Signature and Novel HIF-Regulated Genes.

Molecular Cancer Research. 2018.

GEO Accession: GSE111653
