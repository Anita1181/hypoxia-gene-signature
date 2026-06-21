import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    roc_curve,
    auc,
    accuracy_score,
    classification_report
)

# ==========================================
# LOAD DATA
# ==========================================

expr_df = pd.read_csv(
    r"D:\Main Folder\Misc\Project1_Hypoxia\GSE111653_GilkesSalmonCounts.csv"
)

deg = pd.read_csv(
    r"D:\Main Folder\Misc\Project1_Hypoxia\results\paired_differential_expression.csv"
)

# ==========================================
# TOP 20 GENES
# ==========================================

TOP_N = 20

top_genes = deg.head(TOP_N)["Gene"].tolist()

filtered = expr_df[
    expr_df.iloc[:, 0].isin(top_genes)
]

X = np.log2(filtered.iloc[:, 1:] + 1).T

# ==========================================
# LABELS
# ==========================================

y = np.array([
    1 if "_1percent" in sample else 0
    for sample in X.index
])

# ==========================================
# TRAIN / TEST
# ==========================================

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42,
    stratify=y
)

# ==========================================
# RANDOM FOREST
# ==========================================

rf = RandomForestClassifier(
    n_estimators=500,
    random_state=42
)

rf.fit(X_train, y_train)

# ==========================================
# PREDICTIONS
# ==========================================

pred = rf.predict(X_test)

prob = rf.predict_proba(X_test)[:, 1]

print("\nAccuracy:")
print(accuracy_score(y_test, pred))

print("\nClassification Report:")
print(classification_report(y_test, pred))

# ==========================================
# ROC CURVE
# ==========================================

fpr, tpr, thresholds = roc_curve(
    y_test,
    prob
)

roc_auc = auc(
    fpr,
    tpr
)

print("\nROC-AUC:")
print(round(roc_auc, 4))

# ==========================================
# PLOT
# ==========================================

plt.figure(figsize=(8, 6))

plt.plot(
    fpr,
    tpr,
    linewidth=2,
    label=f"AUC = {roc_auc:.3f}"
)

plt.plot(
    [0, 1],
    [0, 1],
    linestyle="--"
)

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")

plt.title("ROC Curve - Top 20 Hypoxia Genes")

plt.legend()

plt.tight_layout()

plt.savefig(
    r"D:\Main Folder\Misc\Project1_Hypoxia\figures\roc_auc_top20.png",
    dpi=300
)

plt.show()