import pandas as pd

# ----------------------------
# Load data
# ----------------------------
counts = pd.read_csv("counts.csv")

gene_col = counts.columns[0]
sample_cols = counts.columns[1:]

meta = pd.read_csv("metadata_clean.csv")

# Normalize metadata sample column name
if "sample" not in meta.columns and "sample_id" in meta.columns:
    meta = meta.rename(columns={"sample_id": "sample"})

if "sample" not in meta.columns:
    raise ValueError("metadata_clean.csv must contain a 'sample' or 'sample_id' column")

# ----------------------------
# QC metrics per sample
# ----------------------------
library_size = counts[sample_cols].sum(axis=0)
detected_genes = (counts[sample_cols] > 0).sum(axis=0)

total_genes = counts.shape[0]
detection_rate = detected_genes / total_genes

low_count_genes = (counts[sample_cols] <= 1).sum(axis=0)
low_count_rate = low_count_genes / total_genes

qc = pd.DataFrame({
    "sample": sample_cols,
    "library_size": library_size.values,
    "detected_genes": detected_genes.values,
    "detection_rate": detection_rate.values,
    "low_count_rate": low_count_rate.values
})

# Add metadata columns if available
qc = qc.merge(meta, on="sample", how="left")

# Warn if metadata didn't match (common)
if "condition" in qc.columns and qc["condition"].isna().any():
    missing = qc.loc[qc["condition"].isna(), "sample"].tolist()
    print("\nWARNING: These samples from counts.csv are missing in metadata_clean.csv (naming mismatch likely):")
    print(missing)

# Clean types for nicer output
qc["library_size"] = qc["library_size"].astype(int)
qc["detected_genes"] = qc["detected_genes"].astype(int)

# ----------------------------
# Flagging rules (simple but useful)
# ----------------------------
lib_med = qc["library_size"].median()
det_med = qc["detected_genes"].median()

qc["flag_low_library"] = qc["library_size"] < (0.5 * lib_med)
qc["flag_low_detect"] = qc["detected_genes"] < (0.7 * det_med)
qc["flag"] = qc["flag_low_library"] | qc["flag_low_detect"]

# Save outputs
qc.to_csv("qc_metrics_per_sample.csv", index=False)

n_flagged = int(qc["flag"].sum())
flagged_samples = qc.loc[qc["flag"], "sample"].tolist()

with open("qc_metrics_summary.txt", "w", encoding="utf-8") as f:
    f.write("QC SUMMARY\n")
    f.write(f"Total samples: {len(qc)}\n")
    f.write(f"Total genes: {total_genes}\n\n")
    f.write(f"Median library_size: {lib_med}\n")
    f.write(f"Median detected_genes: {det_med}\n\n")
    f.write(f"Flagged samples (n={n_flagged}):\n")
    for s in flagged_samples:
        f.write(f"- {s}\n")

print("Saved: qc_metrics_per_sample.csv")
print("Saved: qc_metrics_summary.txt")
print("\nFlagged samples:", flagged_samples if flagged_samples else "None")
print("\nDAY 9 COMPLETE")
