#imports + file path
from pathlib import Path

fastq_path = Path("sample.fastq")
print("Looking for:", fastq_path.resolve())
print("Exists:", fastq_path.exists())

#FASTQ reader generator

def fastq_iter(path):
    """Yield FASTQ records: (header, seq, plus, qual)."""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
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

#basic stats skeleton
def basic_stats(path, preview_reads=3):
    n = 0
    lengths = []
    preview = []

    for h, s, p, q in fastq_iter(path):
        n += 1
        lengths.append(len(s))

        if n <= preview_reads:
            s50 = s[:50] + ("..." if len(s) > 50 else "")
            preview.append((h, s50, len(s)))

    if n == 0:
        raise ValueError("No reads found. FASTQ empty?")

    avg_len = sum(lengths) / n

    print(f"File: {path}")
    print(f"Reads: {n}")
    print(f"Length min/avg/max: {min(lengths)} / {avg_len:.1f} / {max(lengths)}")
    print("\nPreview:")
    for h, s50, L in preview:
        print(h)
        print("SEQ(50bp):", s50)
        print("LEN:", L)
        print()

#main block
if __name__ == "__main__":
    fastq_path = Path("sample.fastq")
    if not fastq_path.exists():
        print("ERROR: sample.fastq not found in this folder.")
    else:
        basic_stats(fastq_path, preview_reads=3)


