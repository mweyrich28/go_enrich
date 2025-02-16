import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = {
    "Type": ["parallel", "single threaded", "parallel", "single threaded"],
    "Mapping": ["go", "go", "ensembl", "ensembl"],
    "Total Time": [7.869, 24.245, 6.101, 18.077],
    "Reading Input Files": [1.028, 1.074, 0.652, 0.652],
    "Optimizing DAG": [1.563, 1.511, 1.526, 1.618],
    "Enrichment Analysis": [1.956, 11.267, 1.496, 7.984],
    "GO Features": [3.506, 10.309, 2.386, 7.762],
}

df = pd.DataFrame(data)

df_melted = df.melt(id_vars=["Type", "Mapping"], var_name="Stage", value_name="Time")

df_melted["Category"] = df_melted["Type"] + " - " + df_melted["Mapping"]

plt.figure(figsize=(12, 9))
sns.set_palette("colorblind")
sns.set_style("ticks")
ax = sns.barplot(
    x="Stage",
    y="Time",
    hue="Category",
    data=df_melted,
    palette="Set2",
    edgecolor="black",
)
for p in ax.patches:
    if p.get_height() >= 0.2:
        ax.annotate(
            f"{p.get_height():.2f}s",
            (p.get_x() + p.get_width() / 2.0, p.get_height()),
            ha="center",
            va="bottom",
            fontsize=12,
            rotation=90,
            fontweight="bold",
            color="black",
            xytext=(0, 5),
            textcoords="offset points",
        )
plt.xlabel("Processing Stage", fontsize=16)
sns.despine()
plt.ylabel("Time (seconds)", fontsize=16)
plt.title("Comparison of Processing Time by Stage", fontsize=18)
plt.xticks(rotation=45, fontsize=14)
plt.yticks(fontsize=14)
plt.legend(title="Category", title_fontsize=14, fontsize=12)
plt.tight_layout()

plt.savefig("./../report/plots/times.png")
