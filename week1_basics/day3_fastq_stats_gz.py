# Day 3: FASTQ stats + validation (supports .fastq and .fastq.gz)

from pathlib import Path
import gzip


def open_text(path: Path):
    """Open plain text or .gz text file transparently."""
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def fastq_iter(path: Path):
    """Yield FASTQ records: (header, seq, plus, qual)."""
    with open_text(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline()
            p = f.readline()
            q = f.readline()
            if not (s and p and q):
                raise ValueError("FASTQ ended unexpectedly (incomplete record).")
            yield h.rstrip("\n\r"), s.rstrip("\n\r"), p.rstrip("\n\r"), q.rstrip("\n\r")


def validate_record(h: str, s: str, p: str, q: str):
    """Basic FASTQ sanity checks."""
    if not h.startswith("@"):
        raise ValueError(f"Invalid header (does not start with @): {h[:80]}")
    if not p.startswith("+"):
        raise ValueError(f"Invalid plus line (does not start with +): {p[:80]}")
    if len(s) != len(q):
        raise ValueError(f"Seq/Qual length mismatch: seq={len(s)} qual={len(q)} header={h[:80]}")


def basic_stats(path: Path, preview_reads: int = 3, max_reads: int | None = None):
    n = 0
    lengths: list[int] = []
    preview: list[tuple[str, str, int]] = []

    for h, s, p, q in fastq_iter(path):
        validate_record(h, s, p, q)

        n += 1
        L = len(s)
        lengths.append(L)

        if n <= preview_reads:
            s50 = s[:50] + ("..." if L > 50 else "")
            preview.append((h, s50, L))

        if max_reads is not None and n >= max_reads:
            break

    if n == 0:
        raise ValueError("No reads found. FASTQ empty?")

    avg_len = sum(lengths) / n

    print(f"File: {path}")
    print(f"Reads processed: {n}")
    print(f"Read length min/avg/max: {min(lengths)} / {avg_len:.1f} / {max(lengths)}")
    print("\nPreview (first reads):")
    for h, s50, L in preview:
        print(h)
        print("SEQ(50bp):", s50)
        print("LEN:", L)
        print()


if __name__ == "__main__":
    # Update this filename to match what you actually have in week1_basics
    fastq_path = Path("demofile.fastq.gz")

    print("Looking for:", fastq_path.resolve())
    print("Exists:", fastq_path.exists())

    if not fastq_path.exists():
        print("ERROR: FASTQ file not found in this folder.")
        print("Fix: put your .fastq/.fastq.gz here OR change fastq_path in the script.")
    else:
        # Tip: set max_reads=100000 for quick test, then remove to process full file
        basic_stats(fastq_path, preview_reads=3, max_reads=200000)
