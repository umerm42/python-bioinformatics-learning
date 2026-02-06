import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ----------------------------
# Load QC metrics
# ----------------------------
qc = pd.read_csv("qc_metrics_per_sample.csv")

# ----------------------------
# Histogram: Library size
# ----------------------------
plt.figure()
plt.hist(qc["library_size"], bins=10)
plt.xlabel("Library size")
plt.ylabel("Number of samples")
plt.title("Library Size Distribution")
plt.tight_layout()
plt.savefig("qc_library_size_hist.png")
plt.close()

# ----------------------------
# Boxplot: Detection rate
# ----------------------------
plt.figure()
plt.boxplot(qc["detection_rate"], vert=False)
plt.xlabel("Detection rate")
plt.title("Gene Detection Rate")
plt.tight_layout()
plt.savefig("qc_detection_rate_boxplot.png")
plt.close()

# ----------------------------
# PCA on counts
# ----------------------------
counts = pd.read_csv("counts.csv")
gene_col = counts.columns[0]
sample_cols = counts.columns[1:]

X = counts[sample_cols].T  # samples x genes

# Scale before PCA
X_scaled = StandardScaler().fit_transform(X)

pca = PCA(n_components=2)
pcs = pca.fit_transform(X_scaled)

pca_df = pd.DataFrame(
    pcs, columns=["PC1", "PC2"]
)
pca_df["sample"] = sample_cols

pca_df = pca_df.merge(qc[["sample", "condition"]], on="sample", how="left")

# Plot PCA
plt.figure()
for cond in pca_df["condition"].unique():
    sub = pca_df[pca_df["condition"] == cond]
    plt.scatter(sub["PC1"], sub["PC2"], label=cond)

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.title("PCA of Samples (Counts)")
plt.legend()
plt.tight_layout()
plt.savefig("qc_pca_samples.png")
plt.close()

print("Saved plots:")
print("- qc_library_size_hist.png")
print("- qc_detection_rate_boxplot.png")
print("- qc_pca_samples.png")
print("\nDAY 10 COMPLETE")
