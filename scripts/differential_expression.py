import pandas as pd
import numpy as np
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent

FILE_PATH = PROJECT_ROOT / "GSE111653_GilkesSalmonCounts.csv"

df = pd.read_csv(FILE_PATH)

genes = df.iloc[:, 0]
expr = np.log2(df.iloc[:, 1:] + 1)

hypoxia_cols = [c for c in expr.columns if "_1percent" in c]
normoxia_cols = [c for c in expr.columns if "_20percent" in c]

print("Hypoxia samples:", len(hypoxia_cols))
print("Normoxia samples:", len(normoxia_cols))

results = []

for i in range(len(genes)):

    hypoxia_mean = expr.loc[i, hypoxia_cols].mean()
    normoxia_mean = expr.loc[i, normoxia_cols].mean()

    log2_fc = hypoxia_mean - normoxia_mean

    results.append(
        [genes.iloc[i], hypoxia_mean, normoxia_mean, log2_fc]
    )

deg = pd.DataFrame(
    results,
    columns=[
        "Gene",
        "HypoxiaMean",
        "NormoxiaMean",
        "Log2FoldChange"
    ]
)

deg = deg.sort_values(
    "Log2FoldChange",
    ascending=False
)

deg.to_csv(
    PROJECT_ROOT / "results" / "differential_expression.csv",
    index=False
)

print(deg.head(20))
