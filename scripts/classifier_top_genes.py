import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
from pathlib import Path

# -------------------------
# Load expression data
# -------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
expr_df = pd.read_csv(
    PROJECT_ROOT / "GSE111653_GilkesSalmonCounts.csv"
)

# -------------------------
# Load DE results
# -------------------------

deg = pd.read_csv(
    PROJECT_ROOT / "results" / "paired_differential_expression.csv"
)

TOP_N = 200

top_genes = deg.head(TOP_N)["Gene"].tolist()

# -------------------------
# Keep only top genes
# -------------------------

filtered = expr_df[
    expr_df.iloc[:,0].isin(top_genes)
]

X = np.log2(filtered.iloc[:,1:] + 1).T

# -------------------------
# Labels
# -------------------------

y = np.array([
    1 if "_1percent" in sample else 0
    for sample in X.index
])

# -------------------------
# Train/Test Split
# -------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42,
    stratify=y
)

# -------------------------
# Random Forest
# -------------------------

rf = RandomForestClassifier(
    n_estimators=500,
    random_state=42
)

rf.fit(X_train, y_train)

pred = rf.predict(X_test)

print("Accuracy:")
print(accuracy_score(y_test, pred))

print("\nClassification Report:")
print(classification_report(y_test, pred))
