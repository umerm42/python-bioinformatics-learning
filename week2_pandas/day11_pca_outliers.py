import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np

# ----------------------------
# Load inputs
# ----------------------------
qc = pd.read_csv("qc_metrics_per_sample.csv")
counts = pd.read_csv("counts.csv")

gene_col = counts.columns[0]
sample_cols = counts.columns[1:]

# Ensure metadata merge works
if "sample" not in qc.columns:
    raise ValueError("qc_metrics_per_sample.csv must have 'sample' column")

# ----------------------------
# PCA
# ----------------------------
X = counts[sample_cols].T  # samples x genes
X_scaled = StandardScaler().fit_transform(X)

pca = PCA(n_components=2)
pcs = pca.fit_transform(X_scaled)

pca_df = pd.DataFrame(pcs, columns=["PC1", "PC2"])
pca_df["sample"] = sample_cols
pca_df = pca_df.merge(qc[["sample", "condition", "library_size", "detection_rate"]], on="sample", how="left")

# Save explained variance
var1 = pca.explained_variance_ratio_[0] * 100
var2 = pca.explained_variance_ratio_[1] * 100

with open("pca_explained_variance.txt", "w") as f:
    f.write(f"PC1: {var1:.2f}%\n")
    f.write(f"PC2: {var2:.2f}%\n")

# ----------------------------
# Outlier detection in PC space (robust-ish)
# rule: distance from centroid > mean + 2*std
# ----------------------------
centroid = pca_df[["PC1", "PC2"]].mean().values
dist = np.sqrt(((pca_df[["PC1","PC2"]].values - centroid) ** 2).sum(axis=1))
pca_df["pc_distance"] = dist

thr = dist.mean() + 2 * dist.std()
pca_df["flag_outlier"] = pca_df["pc_distance"] > thr

outliers = pca_df[pca_df["flag_outlier"]].copy()
outliers.to_csv("qc_outliers.csv", index=False)

# ----------------------------
# Plot annotated PCA
# ----------------------------
plt.figure()

for cond in pca_df["condition"].dropna().unique():
    sub = pca_df[pca_df["condition"] == cond]
    plt.scatter(sub["PC1"], sub["PC2"], label=str(cond))

# If any samples have missing condition, plot them too
missing = pca_df[pca_df["condition"].isna()]
if len(missing) > 0:
    plt.scatter(missing["PC1"], missing["PC2"], label="missing_condition")

# annotate sample names
for _, row in pca_df.iterrows():
    plt.text(row["PC1"], row["PC2"], str(row["sample"]), fontsize=8)

# highlight outliers
if len(outliers) > 0:
    plt.scatter(outliers["PC1"], outliers["PC2"], marker="x", s=80)

plt.xlabel(f"PC1 ({var1:.1f}%)")
plt.ylabel(f"PC2 ({var2:.1f}%)")
plt.title("Annotated PCA + Outlier Flags")
plt.legend()
plt.tight_layout()
plt.savefig("qc_pca_annotated.png")
plt.close()

print("Saved: qc_pca_annotated.png")
print("Saved: qc_outliers.csv")
print("Saved: pca_explained_variance.txt")
print("Outliers flagged:", outliers["sample"].tolist() if len(outliers) else "None")
print("\nDAY 11 COMPLETE")
