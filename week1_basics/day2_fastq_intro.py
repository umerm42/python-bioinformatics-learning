# Day 2: FASTQ basic stats (supports .fastq and .fastq.gz)

from pathlib import Path
import gzip


def open_text(path: Path):
    """Open plain text or .gz text file transparently."""
    path_str = str(path)
    if path_str.endswith(".gz"):
        return gzip.open(path_str, "rt", encoding="utf-8", errors="replace")
    return open(path_str, "r", encoding="utf-8", errors="replace")


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
            yield h.strip(), s.strip(), p.strip(), q.strip()


def basic_stats(path: Path, preview_reads: int = 3, max_reads: int | None = None):
    n = 0
    lengths: list[int] = []
    preview: list[tuple[str, str, int]] = []

    for h, s, p, q in fastq_iter(path):
        n += 1
        lengths.append(len(s))

        if n <= preview_reads:
            s50 = s[:50] + ("..." if len(s) > 50 else "")
            preview.append((h, s50, len(s)))

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
    # Change this if your filename is different:
    fastq_path = Path("demofile.fastq.gz")

    print("Looking for:", fastq_path.resolve())
    print("Exists:", fastq_path.exists())

    if not fastq_path.exists():
        print("ERROR: demofile.fastq.gz not found in this folder.")
        print("Fix: place the FASTQ in this folder or update fastq_path in the script.")
    else:
        basic_stats(fastq_path, preview_reads=3)
