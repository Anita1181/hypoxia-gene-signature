import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


# ==========================
# LOAD DATA
# ==========================

FILE_PATH = PROJECT_ROOT / "GSE111653_GilkesSalmonCounts.csv"

df = pd.read_csv(FILE_PATH)

print("Original Shape:", df.shape)


# ==========================
# STEP 1
# SEPARATING GENES + EXPRESSION
# ==========================

genes = df.iloc[:, 0]

expression = df.iloc[:, 1:]

print("Expression Matrix Shape:", expression.shape)


# ==========================
# STEP 2
# LOG2 TRANSFORMATION
# ==========================

expression_log = np.log2(expression + 1)

print("Log transform complete")


# ==========================
# STEP 3
# REMOVE LOW EXPRESSION GENES
# ==========================

mean_expression = expression_log.mean(axis=1)

mask = mean_expression > 1

expression_filtered = expression_log[mask]

genes_filtered = genes[mask]

print("Genes before filtering:", len(genes))
print("Genes after filtering :", len(genes_filtered))


# ==========================
# STEP 4
# PCA
# ==========================

X = expression_filtered.T

scaler = StandardScaler()

X_scaled = scaler.fit_transform(X)

pca = PCA(n_components=2)

pca_result = pca.fit_transform(X_scaled)

print("\nPCA Complete")
print("PC1 variance:", round(pca.explained_variance_ratio_[0] * 100, 2), "%")
print("PC2 variance:", round(pca.explained_variance_ratio_[1] * 100, 2), "%")


# ==========================
# STEP 5
# VISUALIZATION
# ==========================

sample_names = expression_filtered.columns

colors = []

for sample in sample_names:
    if "_1percent" in sample:
        colors.append("red")      # Hypoxia
    else:
        colors.append("blue")     # Normoxia

plt.figure(figsize=(12, 8))

plt.scatter(
    pca_result[:, 0],
    pca_result[:, 1],
    c=colors,
    alpha=0.8,
    s=60
)

for i, sample in enumerate(sample_names):
    plt.annotate(
        sample,
        (pca_result[i, 0], pca_result[i, 1]),
        fontsize=7
    )

plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)"
)

plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)"
)

plt.title("PCA of GSE111653 Hypoxia Dataset")

# Custom legend
from matplotlib.lines import Line2D

legend_elements = [
    Line2D(
        [0], [0],
        marker='o',
        color='w',
        label='Hypoxia (1% O₂)',
        markerfacecolor='red',
        markersize=8
    ),
    Line2D(
        [0], [0],
        marker='o',
        color='w',
        label='Normoxia (20% O₂)',
        markerfacecolor='blue',
        markersize=8
    )
]

plt.legend(handles=legend_elements)

plt.tight_layout()

plt.savefig(
    PROJECT_ROOT / "figures" / "pca_hypoxia_vs_normoxia.png",
    dpi=300
)

plt.show()
