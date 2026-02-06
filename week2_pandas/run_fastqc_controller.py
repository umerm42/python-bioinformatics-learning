from pathlib import Path
import subprocess
from datetime import datetime

# SAFE default for shared lab server:
THREADS = 2

WORKDIR = Path(".").resolve()
SAMPLE_FILE = WORKDIR / "sample_run1.txt"
QC_DIR = WORKDIR / "run" / "01_seq_qc"
SUMMARY_TSV = QC_DIR / "fastqc_run_summary.tsv"

def log(msg: str) -> None:
    ts = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    print(ts, msg, flush=True)

def run(cmd) -> None:
    log("RUN: " + " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)

def parse_sample_file(path: Path):
    # tab-separated: sample_id \t r1 \t r2
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 3:
                raise ValueError(
                    f"Line {line_no}: expected 3 TAB-separated columns: sample, r1, r2. Got: {line}"
                )
            sample, r1, r2 = parts
            yield sample, Path(r1), Path(r2)

def main():
    if not SAMPLE_FILE.exists():
        raise FileNotFoundError(f"Not found: {SAMPLE_FILE}")

    QC_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    failures = 0

    log(f"WORKDIR: {WORKDIR}")
    log(f"THREADS: {THREADS}")
    log("Starting FastQC batch...")

    for sample, r1, r2 in parse_sample_file(SAMPLE_FILE):
        outdir = QC_DIR / sample
        outdir.mkdir(parents=True, exist_ok=True)

        # Inputs are filenames relative to WORKDIR (same folder)
        r1_path = WORKDIR / r1
        r2_path = WORKDIR / r2

        if not r1_path.exists():
            reason = f"Missing R1: {r1_path}"
            log(f"[FAIL] {sample}: {reason}")
            rows.append((sample, str(r1_path), str(r2_path), "FAIL", reason))
            failures += 1
            continue

        if not r2_path.exists():
            reason = f"Missing R2: {r2_path}"
            log(f"[FAIL] {sample}: {reason}")
            rows.append((sample, str(r1_path), str(r2_path), "FAIL", reason))
            failures += 1
            continue

        try:
            log(f"FastQC: {sample}")
            run([
                "fastqc",
                str(r1_path), str(r2_path),
                "-t", str(THREADS),
                "-o", str(outdir)
            ])
            rows.append((sample, str(r1_path), str(r2_path), "OK", ""))
        except subprocess.CalledProcessError as e:
            reason = f"fastqc exit {e.returncode}"
            log(f"[FAIL] {sample}: {reason}")
            rows.append((sample, str(r1_path), str(r2_path), "FAIL", reason))
            failures += 1

    # Write summary TSV
    with open(SUMMARY_TSV, "w", encoding="utf-8") as f:
        f.write("sample\tr1\tr2\tstatus\tnote\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    log(f"Saved summary: {SUMMARY_TSV}")

    # MultiQC once (aggregates all FastQC outputs)
    log("Running MultiQC...")
    run(["multiqc", str(QC_DIR), "--force", "-o", str(QC_DIR), "-n", "01_report.html"])

    if failures:
        log(f"DONE with failures: {failures}. See {SUMMARY_TSV}")
    else:
        log("DONE: all samples processed successfully")

if __name__ == "__main__":
    main()
