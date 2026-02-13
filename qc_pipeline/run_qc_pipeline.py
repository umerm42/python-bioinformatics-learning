#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ----------------------------
# Data structures
# ----------------------------
@dataclass
class Sample:
    name: str
    fq1: Path
    fq2: Path


@dataclass
class StepResult:
    sample: str
    step: str
    status: str  # "OK" | "SKIP" | "FAIL"
    detail: str


# ----------------------------
# Utilities
# ----------------------------
def die(msg: str, code: int = 1) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)


def now_iso() -> str:
    return dt.datetime.now().replace(microsecond=0).isoformat()


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def file_ok(p: Path, min_bytes: int = 1) -> bool:
    return p.exists() and p.is_file() and p.stat().st_size >= min_bytes


def run_cmd(
    cmd: List[str],
    log_path: Path,
    dry_run: bool = False,
    cwd: Optional[Path] = None,
) -> None:
    """
    Run command, append stdout/stderr to log.
    Fail fast if command returns non-zero.
    """
    cmd_str = " ".join([shlex_quote(x) for x in cmd])
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n[{now_iso()}] $ {cmd_str}\n")
        log.flush()

        if dry_run:
            log.write("[DRY-RUN] Command not executed.\n")
            return

        proc = subprocess.run(
            cmd,
            stdout=log,
            stderr=log,
            cwd=str(cwd) if cwd else None,
            check=False,
            text=True,
        )
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)


def shlex_quote(s: str) -> str:
    # minimal safe quoting for logs
    if not s:
        return "''"
    if any(ch.isspace() or ch in "\"'()[]{}$&;|<>*" for ch in s):
        return "'" + s.replace("'", "'\"'\"'") + "'"
    return s


def which_or_die(tool_name: str, tool_cmd: str) -> str:
    """
    Resolve executable. tool_cmd may already be a path or a command in PATH.
    """
    resolved = shutil.which(tool_cmd)
    if not resolved:
        die(
            f"Tool '{tool_name}' not found: '{tool_cmd}'. "
            f"Fix PATH or config.yaml tools.{tool_name}"
        )
    return resolved


def get_tool_version(cmd: str, log_path: Path) -> str:
    """
    Best-effort version capture. Not all tools behave consistently.
    """
    candidates = [
        [cmd, "--version"],
        [cmd, "-v"],
        [cmd, "version"],
    ]
    for c in candidates:
        try:
            out = subprocess.run(
                c, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False
            ).stdout.strip()
            if out:
                with log_path.open("a", encoding="utf-8") as log:
                    log.write(f"\n[{now_iso()}] VERSION {cmd}: {' '.join(c[1:])}\n{out}\n")
                return out.splitlines()[0][:200]
        except Exception:
            continue
    return "unknown"


def load_yaml_config(path: Path) -> Dict:
    """
    Requires PyYAML. That’s normal for real work.
    """
    try:
        import yaml  # type: ignore
    except ImportError:
        die(
            "Missing dependency: PyYAML.\n"
            "Install it once: pip install pyyaml\n"
            "Or in conda: conda install pyyaml"
        )
    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        die("config.yaml did not parse into a dictionary.")
    return cfg


def load_samples_tsv(path: Path) -> List[Sample]:
    if not path.exists():
        die(f"Samples file not found: {path}")
    samples: List[Sample] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"sample", "fq1", "fq2"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            die(f"samples.tsv must have columns: sample, fq1, fq2 (tab-separated). Got: {reader.fieldnames}")
        for row in reader:
            name = (row.get("sample") or "").strip()
            fq1 = (row.get("fq1") or "").strip()
            fq2 = (row.get("fq2") or "").strip()
            if not name or not fq1 or not fq2:
                die(f"Bad row in samples.tsv (missing fields): {row}")
            samples.append(Sample(name=name, fq1=Path(fq1), fq2=Path(fq2)))
    # uniqueness
    seen = set()
    for s in samples:
        if s.name in seen:
            die(f"Duplicate sample name in samples.tsv: {s.name}")
        seen.add(s.name)
    return samples


def validate_inputs(samples: List[Sample]) -> None:
    for s in samples:
        if not s.fq1.exists():
            die(f"Missing FASTQ for sample {s.name}: {s.fq1}")
        if not s.fq2.exists():
            die(f"Missing FASTQ for sample {s.name}: {s.fq2}")


# ----------------------------
# Pipeline steps
# ----------------------------
def step_fastqc(
    fastqc_cmd: str,
    fq1: Path,
    fq2: Path,
    outdir: Path,
    threads: int,
    log_path: Path,
    dry_run: bool,
) -> Tuple[Path, Path]:
    """
    Run FastQC on a pair of reads.
    Returns expected .zip outputs.
    """
    ensure_dir(outdir)

    # FastQC naming convention uses input filename stem + _fastqc.zip
    zip1 = outdir / f"{fq1.name}_fastqc.zip"
    zip2 = outdir / f"{fq2.name}_fastqc.zip"

    if file_ok(zip1) and file_ok(zip2):
        return zip1, zip2

    cmd = [
        fastqc_cmd,
        "--threads", str(threads),
        "--outdir", str(outdir),
        str(fq1),
        str(fq2),
    ]
    run_cmd(cmd, log_path=log_path, dry_run=dry_run)
    return zip1, zip2


def step_fastp(
    fastp_cmd: str,
    sample: str,
    fq1: Path,
    fq2: Path,
    outdir: Path,
    threads: int,
    extra_args: str,
    log_path: Path,
    dry_run: bool,
) -> Tuple[Path, Path, Path, Path]:
    """
    Run fastp on paired-end reads.
    Returns (trim1, trim2, json, html).
    """
    sdir = outdir / sample
    ensure_dir(sdir)

    trim1 = sdir / f"{sample}_R1.trim.fq.gz"
    trim2 = sdir / f"{sample}_R2.trim.fq.gz"
    jsonp = sdir / f"{sample}.fastp.json"
    htmlp = sdir / f"{sample}.fastp.html"

    if file_ok(trim1) and file_ok(trim2) and file_ok(jsonp) and file_ok(htmlp):
        return trim1, trim2, jsonp, htmlp

    cmd = [
        fastp_cmd,
        "--in1", str(fq1),
        "--in2", str(fq2),
        "--out1", str(trim1),
        "--out2", str(trim2),
        "--json", str(jsonp),
        "--html", str(htmlp),
        "--thread", str(threads),
    ]

    # split extra args safely
    if extra_args.strip():
        cmd.extend(extra_args.strip().split())

    run_cmd(cmd, log_path=log_path, dry_run=dry_run)
    return trim1, trim2, jsonp, htmlp


def step_multiqc(
    multiqc_cmd: str,
    scan_dir: Path,
    outdir: Path,
    log_path: Path,
    dry_run: bool,
) -> Path:
    ensure_dir(outdir)
    report = outdir / "multiqc_report.html"
    if file_ok(report):
        return report

    cmd = [
        multiqc_cmd,
        str(scan_dir),
        "--outdir", str(outdir),
        "--force",
    ]
    run_cmd(cmd, log_path=log_path, dry_run=dry_run)
    return report


def write_report(
    outdir: Path,
    cfg: Dict,
    samples: List[Sample],
    results: List[StepResult],
    tool_versions: Dict[str, str],
) -> Path:
    report_path = outdir / "QC_REPORT.md"

    total = len(samples)
    failed = [r for r in results if r.status == "FAIL"]
    ok = [r for r in results if r.status == "OK"]
    skipped = [r for r in results if r.status == "SKIP"]

    lines: List[str] = []
    lines.append("# RNA-seq QC Report\n")
    lines.append(f"- Generated: {now_iso()}\n")
    lines.append(f"- Samples: {total}\n")
    lines.append(f"- Outdir: `{outdir}`\n")

    lines.append("\n## Tool Versions\n")
    for k, v in tool_versions.items():
        lines.append(f"- **{k}**: {v}\n")

    lines.append("\n## Steps Enabled\n")
    steps = cfg.get("steps", {})
    for k in ["fastqc_raw", "trim_fastp", "fastqc_trimmed", "multiqc"]:
        lines.append(f"- {k}: {bool(steps.get(k, False))}\n")

    lines.append("\n## Key Outputs\n")
    lines.append(f"- Raw FastQC: `{outdir / '01_fastqc_raw'}`\n")
    lines.append(f"- fastp trimmed reads: `{outdir / '02_trim_fastp'}`\n")
    lines.append(f"- Trimmed FastQC: `{outdir / '03_fastqc_trimmed'}`\n")
    lines.append(f"- MultiQC: `{outdir / '04_multiqc' / 'multiqc_report.html'}`\n")

    lines.append("\n## Run Summary\n")
    lines.append(f"- OK: {len(ok)}\n")
    lines.append(f"- Skipped: {len(skipped)}\n")
    lines.append(f"- Failed: {len(failed)}\n")

    if failed:
        lines.append("\n## Failures (must fix)\n")
        for r in failed:
            lines.append(f"- **{r.sample}** / {r.step}: {r.detail}\n")

    lines.append("\n## Interpretation Notes\n")
    lines.append("- Open `04_multiqc/multiqc_report.html` and scan for outliers.\n")
    lines.append("- If trimming removes a large fraction of reads, investigate adapter/quality issues.\n")
    lines.append("- If one sample is a severe outlier across metrics, confirm sample identity / contamination.\n")

    report_path.write_text("".join(lines), encoding="utf-8")
    return report_path


# ----------------------------
# Main
# ----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="End-to-end RNA-seq QC pipeline (FastQC + fastp + MultiQC) with report generation."
    )
    p.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    p.add_argument("--limit", type=int, default=0, help="Run only first N samples (0 = all)")
    p.add_argument("--dry-run", action="store_true", help="Print/log commands but do not execute")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        die(f"Config file not found: {cfg_path}")

    cfg = load_yaml_config(cfg_path)

    samples_tsv = Path(cfg.get("samples_tsv", "samples.tsv"))
    outdir = Path(cfg.get("outdir", "results"))
    threads = int(cfg.get("threads", 4))

    tools = cfg.get("tools", {}) or {}
    steps = cfg.get("steps", {}) or {}

    ensure_dir(outdir)
    ensure_dir(outdir / "logs")
    pipeline_log = outdir / "logs" / "pipeline.log"

    # Resolve tools early (fail fast)
    fastqc_cmd = which_or_die("fastqc", str(tools.get("fastqc", "fastqc")))
    fastp_cmd = which_or_die("fastp", str(tools.get("fastp", "fastp")))
    multiqc_cmd = which_or_die("multiqc", str(tools.get("multiqc", "multiqc")))

    # Capture versions (best-effort)
    tool_versions = {
        "fastqc": get_tool_version(fastqc_cmd, pipeline_log),
        "fastp": get_tool_version(fastp_cmd, pipeline_log),
        "multiqc": get_tool_version(multiqc_cmd, pipeline_log),
        "python": sys.version.split()[0],
    }

    print("\n=== RUNNING QC PIPELINE ===")
    print(f"Config: {cfg_path}")
    print(f"Samples: {samples_tsv}")
    print(f"Outdir: {outdir}")
    print(f"Threads: {threads}")
    if args.dry_run:
        print("Mode: DRY-RUN (no commands executed)")

    samples = load_samples_tsv(samples_tsv)
    if args.limit and args.limit > 0:
        samples = samples[: args.limit]

    validate_inputs(samples)

    results: List[StepResult] = []

    # Step output roots
    d_fastqc_raw = outdir / "01_fastqc_raw"
    d_trim = outdir / "02_trim_fastp"
    d_fastqc_trim = outdir / "03_fastqc_trimmed"
    d_multiqc = outdir / "04_multiqc"

    fastp_extra = (cfg.get("fastp", {}) or {}).get("extra_args", "")

    for s in samples:
        print(f"\n--- Sample: {s.name} ---")

        # FastQC raw
        if steps.get("fastqc_raw", True):
            try:
                zip1 = d_fastqc_raw / f"{s.fq1.name}_fastqc.zip"
                zip2 = d_fastqc_raw / f"{s.fq2.name}_fastqc.zip"
                if file_ok(zip1) and file_ok(zip2):
                    results.append(StepResult(s.name, "fastqc_raw", "SKIP", "outputs exist"))
                else:
                    step_fastqc(
                        fastqc_cmd=fastqc_cmd,
                        fq1=s.fq1,
                        fq2=s.fq2,
                        outdir=d_fastqc_raw,
                        threads=threads,
                        log_path=pipeline_log,
                        dry_run=args.dry_run,
                    )
                    results.append(StepResult(s.name, "fastqc_raw", "OK", ""))
            except Exception as e:
                results.append(StepResult(s.name, "fastqc_raw", "FAIL", str(e)))
                continue  # don't proceed for this sample

        # fastp trim
        trim1 = trim2 = None
        if steps.get("trim_fastp", True):
            try:
                sdir = d_trim / s.name
                trim1p = sdir / f"{s.name}_R1.trim.fq.gz"
                trim2p = sdir / f"{s.name}_R2.trim.fq.gz"
                jsonp = sdir / f"{s.name}.fastp.json"
                htmlp = sdir / f"{s.name}.fastp.html"
                if file_ok(trim1p) and file_ok(trim2p) and file_ok(jsonp) and file_ok(htmlp):
                    results.append(StepResult(s.name, "trim_fastp", "SKIP", "outputs exist"))
                    trim1, trim2 = trim1p, trim2p
                else:
                    trim1, trim2, _, _ = step_fastp(
                        fastp_cmd=fastp_cmd,
                        sample=s.name,
                        fq1=s.fq1,
                        fq2=s.fq2,
                        outdir=d_trim,
                        threads=threads,
                        extra_args=str(fastp_extra),
                        log_path=pipeline_log,
                        dry_run=args.dry_run,
                    )
                    results.append(StepResult(s.name, "trim_fastp", "OK", ""))
            except Exception as e:
                results.append(StepResult(s.name, "trim_fastp", "FAIL", str(e)))
                continue

        # FastQC trimmed
        if steps.get("fastqc_trimmed", True):
            if trim1 is None or trim2 is None:
                # if trimming disabled, you can choose to fastqc raw only; we won't guess.
                results.append(StepResult(s.name, "fastqc_trimmed", "SKIP", "no trimmed reads"))
            else:
                try:
                    zip1 = d_fastqc_trim / f"{Path(trim1).name}_fastqc.zip"
                    zip2 = d_fastqc_trim / f"{Path(trim2).name}_fastqc.zip"
                    if file_ok(zip1) and file_ok(zip2):
                        results.append(StepResult(s.name, "fastqc_trimmed", "SKIP", "outputs exist"))
                    else:
                        step_fastqc(
                            fastqc_cmd=fastqc_cmd,
                            fq1=Path(trim1),
                            fq2=Path(trim2),
                            outdir=d_fastqc_trim,
                            threads=threads,
                            log_path=pipeline_log,
                            dry_run=args.dry_run,
                        )
                        results.append(StepResult(s.name, "fastqc_trimmed", "OK", ""))
                except Exception as e:
                    results.append(StepResult(s.name, "fastqc_trimmed", "FAIL", str(e)))
                    continue

    # MultiQC (once)
    if steps.get("multiqc", True):
        try:
            report = d_multiqc / "multiqc_report.html"
            if file_ok(report):
                results.append(StepResult("ALL", "multiqc", "SKIP", "report exists"))
            else:
                step_multiqc(
                    multiqc_cmd=multiqc_cmd,
                    scan_dir=outdir,
                    outdir=d_multiqc,
                    log_path=pipeline_log,
                    dry_run=args.dry_run,
                )
                results.append(StepResult("ALL", "multiqc", "OK", ""))
        except Exception as e:
            results.append(StepResult("ALL", "multiqc", "FAIL", str(e)))

    report_path = write_report(outdir=outdir, cfg=cfg, samples=samples, results=results, tool_versions=tool_versions)

    # Print final status
    failed = [r for r in results if r.status == "FAIL"]
    print("\n=== DONE ===")
    print(f"Report: {report_path}")
    print(f"Log: {pipeline_log}")
    if failed:
        print("\nFailures:")
        for r in failed:
            print(f"- {r.sample} / {r.step}: {r.detail}")
        sys.exit(2)
    else:
        print("Pipeline executed successfully.")


if __name__ == "__main__":
    main()
