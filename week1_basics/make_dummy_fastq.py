import random

def rand_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))

def main():
    out = "sample.fastq"
    n_reads = 20
    read_len = 75

    with open(out, "w", encoding="utf-8") as f:
        for i in range(n_reads):
            seq = rand_seq(read_len)
            qual = "I" * read_len  # dummy high-quality Phred33
            f.write(f"@read{i}\n{seq}\n+\n{qual}\n")

    print(f"Created {out} with {n_reads} reads of length {read_len}")

if __name__ == "__main__":
    main()
