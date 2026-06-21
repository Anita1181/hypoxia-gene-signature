import pandas as pd
import numpy as np

FILE_PATH = r"D:\Main Folder\Misc\Project1_Hypoxia\GSE111653_GilkesSalmonCounts.csv"

df = pd.read_csv(FILE_PATH)

genes = df.iloc[:, 0]

expr = np.log2(df.iloc[:, 1:] + 1)

# -------------------------
# find cell lines
# -------------------------

cell_lines = set()

for col in expr.columns:
    if "_1percent" in col:
        cell_lines.add(col.replace("_1percent", ""))

cell_lines = sorted(cell_lines)

print("Matched cell lines:", len(cell_lines))

# -------------------------
# calculate paired fold changes
# -------------------------

fold_changes = []

for cell_line in cell_lines:

    hypoxia_col = f"{cell_line}_1percent"
    normoxia_col = f"{cell_line}_20percent"

    fc = expr[hypoxia_col] - expr[normoxia_col]

    fold_changes.append(fc)

fold_changes = pd.concat(fold_changes, axis=1)

fold_changes.columns = cell_lines

# -------------------------
# Average response
# -------------------------

mean_fc = fold_changes.mean(axis=1)

std_fc = fold_changes.std(axis=1)

results = pd.DataFrame({
    "Gene": genes,
    "MeanLog2FC": mean_fc,
    "StdLog2FC": std_fc
})

results = results.sort_values(
    by="MeanLog2FC",
    ascending=False
)

results.to_csv(
    r"D:\Main Folder\Misc\Project1_Hypoxia\results\paired_differential_expression.csv",
    index=False
)

print("\nTop 20 Upregulated Genes\n")

print(
    results[
        ["Gene", "MeanLog2FC"]
    ].head(20)
)

print("\nSaved to results folder.")