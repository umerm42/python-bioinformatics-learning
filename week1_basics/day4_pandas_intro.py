import pandas as pd

# Load table
df = pd.read_csv("counts.csv")

print("\n--- HEAD ---")
print(df.head())

print("\n--- INFO ---")
print(df.info())

# Basic stats
print("\n--- DESCRIBE ---")
print(df.describe())

# Filter: keep genes with total counts > 50
df["sum_counts"] = df.iloc[:, 1:].sum(axis=1)
filtered = df[df["sum_counts"] > 50]

print("\nFiltered genes:", filtered.shape[0])

# Save filtered table
filtered.to_csv("counts_filtered.csv", index=False)

print("\nSaved: counts_filtered.csv")
