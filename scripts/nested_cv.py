import pandas as pd
import numpy as np
from pathlib import Path

from sklearn.model_selection import StratifiedGroupKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score

# -------------------------
# Load expression data
# -------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
expr_df = pd.read_csv(PROJECT_ROOT / "GSE111653_GilkesSalmonCounts.csv")

genes = expr_df.iloc[:, 0]
expr = expr_df.iloc[:, 1:]

sample_names = expr.columns.tolist()

labels = np.array([1 if "_1percent" in s else 0 for s in sample_names])
groups = np.array([s.replace("_1percent", "").replace("_20percent", "") for s in sample_names])

X_full = expr.T
X_full.index = sample_names

TOP_N = 200     
N_FOLDS = 5

sgkf = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)

fold_results = []

print(f"\nRunning Nested {N_FOLDS}-Fold Cross-Validation (TOP_N = {TOP_N})")
print("=" * 70)

for fold_num, (train_idx, test_idx) in enumerate(
    sgkf.split(X_full, labels, groups), start=1
):
    train_samples = [sample_names[i] for i in train_idx]
    test_samples = [sample_names[i] for i in test_idx]
    train_groups = sorted(set(groups[train_idx]))

    fold_changes = []
    for cell_line in train_groups:
        hyp_col = f"{cell_line}_1percent"
        norm_col = f"{cell_line}_20percent"
        if hyp_col in expr.columns and norm_col in expr.columns:
            log_hyp = np.log2(expr[hyp_col] + 1)
            log_norm = np.log2(expr[norm_col] + 1)
            fold_changes.append(log_hyp - log_norm)

    fold_changes = pd.concat(fold_changes, axis=1)
    mean_fc = fold_changes.mean(axis=1)

    de_results = pd.DataFrame({"Gene": genes, "MeanLog2FC": mean_fc})
    de_results["AbsFC"] = de_results["MeanLog2FC"].abs()
    top_genes_idx = de_results.sort_values("AbsFC", ascending=False).head(TOP_N).index

    X_train = np.log2(X_full.loc[train_samples].iloc[:, top_genes_idx] + 1)
    X_test = np.log2(X_full.loc[test_samples].iloc[:, top_genes_idx] + 1)
    y_train = labels[train_idx]
    y_test = labels[test_idx]

    rf = RandomForestClassifier(n_estimators=500, random_state=42)
    rf.fit(X_train, y_train)

    pred = rf.predict(X_test)
    pred_proba = rf.predict_proba(X_test)[:, 1]

    acc = accuracy_score(y_test, pred)
    auc = roc_auc_score(y_test, pred_proba)

    fold_results.append({
        "fold": fold_num,
        "accuracy": acc,
        "auc": auc,
        "n_test_cell_lines": len(set(groups[test_idx]))
    })

    print(f"Fold {fold_num}  |  Accuracy: {acc:.3f}  |  ROC-AUC: {auc:.3f}  "
          f"|  Test cell lines: {len(set(groups[test_idx]))}")

accs = [r["accuracy"] for r in fold_results]
aucs = [r["auc"] for r in fold_results]

print("=" * 70)
print(f"{'Mean Accuracy:':<20}{np.mean(accs):.3f}  (± {np.std(accs):.3f})")
print(f"{'Mean ROC-AUC:':<20}{np.mean(aucs):.3f}  (± {np.std(aucs):.3f})")
print("=" * 70)
