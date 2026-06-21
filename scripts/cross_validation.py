import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score

expr_df = pd.read_csv(
    r"D:\Main Folder\Misc\Project1_Hypoxia\GSE111653_GilkesSalmonCounts.csv"
)

deg = pd.read_csv(
    r"D:\Main Folder\Misc\Project1_Hypoxia\results\paired_differential_expression.csv"
)

TOP_N = 20

top_genes = deg.head(TOP_N)["Gene"].tolist()

filtered = expr_df[
    expr_df.iloc[:, 0].isin(top_genes)
]

X = np.log2(filtered.iloc[:, 1:] + 1).T

y = np.array([
    1 if "_1percent" in sample else 0
    for sample in X.index
])

rf = RandomForestClassifier(
    n_estimators=500,
    random_state=42
)

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)

accuracy_scores = cross_val_score(
    rf,
    X,
    y,
    cv=cv,
    scoring="accuracy"
)

auc_scores = cross_val_score(
    rf,
    X,
    y,
    cv=cv,
    scoring="roc_auc"
)

print("Accuracy Scores:")
print(accuracy_scores)

print("\nMean Accuracy:")
print(accuracy_scores.mean())

print("\nStd Accuracy:")
print(accuracy_scores.std())

print("\nROC-AUC Scores:")
print(auc_scores)

print("\nMean ROC-AUC:")
print(auc_scores.mean())