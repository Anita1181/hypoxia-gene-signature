import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent

df = pd.read_csv(
    PROJECT_ROOT / "results" / "paired_differential_expression.csv"
)

top20 = df.head(20)

plt.figure(figsize=(10, 8))

plt.barh(
    top20["Gene"],
    top20["MeanLog2FC"]
)

plt.xlabel("Average Log2 Fold Change")

plt.title("Top 20 Hypoxia-Induced Genes")

plt.gca().invert_yaxis()

plt.tight_layout()

plt.savefig(
    PROJECT_ROOT / "figures" / "top20_hypoxia_genes.png",
    dpi=300
)

plt.show()
