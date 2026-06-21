import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    r"D:\Main Folder\Misc\Project1_Hypoxia\results\paired_differential_expression.csv"
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
    r"D:\Main Folder\Misc\Project1_Hypoxia\figures\top20_hypoxia_genes.png",
    dpi=300
)

plt.show()