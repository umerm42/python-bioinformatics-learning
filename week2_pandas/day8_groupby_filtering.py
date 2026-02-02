import pandas as pd

# ----------------------------
# Load RNA-seq count table
# ----------------------------
counts = pd.read_csv("counts.csv")

print("\nCounts shape:", counts.shape)
print(counts.head())

# ----------------------------
# Sample-wise library size
# ----------------------------
# sums each sample column (all columns except gene)
lib_sizes = counts.iloc[:, 1:].sum()
lib_sizes = lib_sizes.reset_index()
lib_sizes.columns = ["sample", "library_size"]

print("\nLibrary sizes:")
print(lib_sizes)

lib_sizes.to_csv("qc_library_sizes.csv", index=False)

# ----------------------------
# Gene-wise filtering
# Rule: expressed (count > 5) in at least 2 samples
# ----------------------------
expr_mask = (counts.iloc[:, 1:] > 5).sum(axis=1) >= 2
filtered = counts.loc[expr_mask].copy()

print("\nGenes before filtering:", counts.shape[0])
print("Genes after filtering:", filtered.shape[0])

filtered.to_csv("counts_expr_filtered.csv", index=False)

# ----------------------------
# Grouped aggregation (mean per condition)
# ----------------------------
meta = pd.read_csv("metadata_clean.csv")

# Your Week 1 metadata likely uses "sample_id"
if "sample" not in meta.columns and "sample_id" in meta.columns:
    meta = meta.rename(columns={"sample_id": "sample"})

required = {"sample", "condition"}
missing = required - set(meta.columns)
if missing:
    raise ValueError(f"metadata_clean.csv missing columns: {missing}. Found: {list(meta.columns)}")

# reshape counts to long format: gene, sample, count
long = filtered.melt(id_vars="gene", var_name="sample", value_name="count")

# merge metadata
merged = long.merge(meta, on="sample", how="left")

# IMPORTANT: detect samples in counts that are not in metadata
n_missing_cond = merged["condition"].isna().sum()
if n_missing_cond > 0:
    # show which samples failed to merge (usually naming mismatch)
    bad_samples = merged.loc[merged["condition"].isna(), "sample"].unique()
    print("\nWARNING: Some samples from counts.csv are missing in metadata_clean.csv")
    print("Missing samples:", list(bad_samples))
    print("These rows will be excluded from group aggregation.\n")

# drop rows without condition for grouping
merged_ok = merged.dropna(subset=["condition"])

grouped = (
    merged_ok
    .groupby(["gene", "condition"], as_index=False)["count"]
    .mean()
)

print("\nGrouped means:")
print(grouped.head())

grouped.to_csv("counts_grouped_mean.csv", index=False)

print("\nDAY 8 PIPELINE COMPLETE")
