import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier

expr_df = pd.read_csv(
    r"D:\Main Folder\Misc\Project1_Hypoxia\GSE111653_GilkesSalmonCounts.csv"
)

deg = pd.read_csv(
    r"D:\Main Folder\Misc\Project1_Hypoxia\results\paired_differential_expression.csv"
)

TOP_N = 20

top_genes = deg.head(TOP_N)["Gene"].tolist()

filtered = expr_df[
    expr_df.iloc[:,0].isin(top_genes)
]

X = np.log2(filtered.iloc[:,1:] + 1).T

y = np.array([
    1 if "_1percent" in sample else 0
    for sample in X.index
])

rf = RandomForestClassifier(
    n_estimators=500,
    random_state=42
)

rf.fit(X, y)

importance_df = pd.DataFrame({
    "Gene": top_genes,
    "Importance": rf.feature_importances_
})

importance_df = importance_df.sort_values(
    by="Importance",
    ascending=False
)

print(importance_df)

importance_df.to_csv(
    r"D:\Main Folder\Misc\Project1_Hypoxia\results\feature_importance.csv",
    index=False
)