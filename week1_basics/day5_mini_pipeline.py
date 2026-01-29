from pathlib import Path
import gzip
import pandas as pd

# ---------------- FASTQ QC ----------------

def open_maybe_gz(path: Path):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")

def fastq_qc(path: Path, max_reads: int = 200_000):
    """
    Returns: (reads_processed, min_len, avg_len, max_len)
    Reads at most max_reads for speed.
    """
    if not path.exists():
        raise FileNotFoundError(f"FASTQ not found: {path}")

    n = 0
    lengths = []

    with open_maybe_gz(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline()
            p = f.readline()
            q = f.readline()
            if not (s and p and q):
                raise ValueError("FASTQ ended unexpectedly (incomplete record).")

            # simple validation
            if not h.startswith("@") or not p.startswith("+"):
                raise ValueError(f"Invalid FASTQ record near header: {h.strip()[:80]}")
            if len(s.strip()) != len(q.strip()):
                raise ValueError(f"Seq/Qual length mismatch near header: {h.strip()[:80]}")

            n += 1
            lengths.append(len(s.strip()))
            if n >= max_reads:
                break

    if n == 0:
        raise ValueError("No reads found (empty FASTQ?)")

    avg_len = sum(lengths) / n
    return n, min(lengths), avg_len, max(lengths)

# ---------------- METADATA CLEAN ----------------

REQUIRED_COLS = ["sample_id", "condition"]

def clean_sample_id(s: str) -> str:
    s = s.strip().replace(" ", "_").replace("-", "_")
    return "".join(c for c in s if c.isalnum() or c == "_")

def clean_metadata(in_csv: Path, out_csv: Path) -> int:
    if not in_csv.exists():
        raise FileNotFoundError(f"Metadata not found: {in_csv}")

    df = pd.read_csv(in_csv)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns {missing}. Found: {list(df.columns)}")

    df["sample_id"] = df["sample_id"].astype(str).apply(clean_sample_id)

    # check empty + duplicates after cleaning
    if (df["sample_id"] == "").any():
        raise ValueError("Found empty sample_id after cleaning.")
    if df["sample_id"].duplicated().any():
        dups = df.loc[df["sample_id"].duplicated(), "sample_id"].tolist()
        raise ValueError(f"Duplicate sample_id after cleaning: {dups}")

    df.to_csv(out_csv, index=False)
    return df.shape[0]

# ---------------- COUNT TABLE FILTER ----------------

def filter_counts(count_csv: Path, out_csv: Path, min_sum: int = 50) -> int:
    if not count_csv.exists():
        raise FileNotFoundError(f"Counts table not found: {count_csv}")

    df = pd.read_csv(count_csv)

    if df.shape[1] < 2:
        raise ValueError("Counts table must have at least 2 columns: gene + >=1 sample")

    # assume first col is gene ID, rest are samples
    df["sum_counts"] = df.iloc[:, 1:].sum(axis=1)
    filtered = df[df["sum_counts"] > min_sum].copy()

    filtered.to_csv(out_csv, index=False)
    return filtered.shape[0]

# ---------------- PIPELINE ----------------

def main():
    print("\n--- MINI RNA-seq PREPROCESSING PIPELINE ---\n")

    # Use small demo file names (recommended)
    fastq = Path("sample.fastq")   # change if needed
    meta_raw = Path("metadata_raw.csv")
    counts = Path("counts.csv")

    # FASTQ QC (optional if file exists)
    if fastq.exists():
        print("Running FASTQ QC...")
        n, mn, avg, mx = fastq_qc(fastq, max_reads=200_000)
        print(f"Reads processed: {n}")
        print(f"Length min/avg/max: {mn}/{avg:.1f}/{mx}\n")
    else:
        print(f"FASTQ step skipped (not found): {fastq}\n")

    # Metadata cleaning
    print("Cleaning metadata...")
    n_meta = clean_metadata(meta_raw, Path("metadata_clean.csv"))
    print(f"Samples: {n_meta}\n")

    # Count filtering
    print("Filtering counts...")
    n_genes = filter_counts(counts, Path("counts_filtered.csv"), min_sum=50)
    print(f"Genes retained: {n_genes}\n")

    print("PIPELINE FINISHED SUCCESSFULLY")

if __name__ == "__main__":
    main()
