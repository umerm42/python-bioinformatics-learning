import csv
from pathlib import Path

REQUIRED_COLS = ["sample_id", "condition"]

def clean_sample_id(s: str) -> str:
    s = s.strip()
    s = s.replace(" ", "_")
    s = s.replace("-", "_")
    # keep it simple: letters, digits, underscore only
    out = []
    for ch in s:
        if ch.isalnum() or ch == "_":
            out.append(ch)
    return "".join(out)

def validate_columns(header):
    missing = [c for c in REQUIRED_COLS if c not in header]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Found: {header}")

def clean_metadata(in_csv: Path, out_csv: Path):
    with open(in_csv, "r", newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames or []
        validate_columns(header)

        rows = []
        seen = set()

        for row in reader:
            sid_raw = row["sample_id"]
            sid = clean_sample_id(sid_raw)

            if not sid:
                raise ValueError(f"Empty sample_id after cleaning. Original: {sid_raw}")

            if sid in seen:
                raise ValueError(f"Duplicate sample_id after cleaning: {sid}")
            seen.add(sid)

            row["sample_id"] = sid
            rows.append(row)

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Cleaned rows: {len(rows)}")
    print(f"Saved: {out_csv}")

if __name__ == "__main__":
    inp = Path("metadata_raw.csv")
    out = Path("metadata_clean.csv")

    if not inp.exists():
        raise FileNotFoundError(
            f"{inp} not found. Create it (see template below) and re-run."
        )

    clean_metadata(inp, out)
