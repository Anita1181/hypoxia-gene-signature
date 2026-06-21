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

## Methods

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
- Top differentially expressed genes (classifier_top_genes.py)

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
- Results should be considered exploratory rather than clinical

---

## References

Ye IC, Fertig EJ, DiGiacomo JW, et al.

Molecular Portrait of Hypoxia in Breast Cancer: A Prognostic Signature and Novel HIF-Regulated Genes.

Molecular Cancer Research. 2018.

GEO Accession: GSE111653