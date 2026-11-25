"""Summarize SLURM array-job progress from sacct/squeue output.

When you launch a Glide virtual-screening array (``glide-array.sh``) you end up
with hundreds of array tasks. This parses the tabular text from ``sacct`` (or
``squeue``) and reports how many tasks are PENDING / RUNNING / COMPLETED /
FAILED, so you can see at a glance whether a big screen finished cleanly.

Schrodinger-free and SLURM-free: it parses a text sample, so it is fully
testable offline. On a real cluster, pipe live output in, e.g.:

    sacct -j <jobid> --format=JobID,State,Elapsed -P | python squeue_summary.py -
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter

# Map raw SLURM state strings to a coarse bucket.
STATE_BUCKETS = {
    "PENDING": "PENDING",
    "PD": "PENDING",
    "RUNNING": "RUNNING",
    "R": "RUNNING",
    "COMPLETING": "RUNNING",
    "CG": "RUNNING",
    "COMPLETED": "COMPLETED",
    "CD": "COMPLETED",
    "FAILED": "FAILED",
    "F": "FAILED",
    "CANCELLED": "FAILED",
    "CA": "FAILED",
    "TIMEOUT": "FAILED",
    "TO": "FAILED",
    "OUT_OF_MEMORY": "FAILED",
    "OOM": "FAILED",
    "NODE_FAIL": "FAILED",
    "NF": "FAILED",
}

BUCKET_ORDER = ["PENDING", "RUNNING", "COMPLETED", "FAILED", "OTHER"]


def _normalize_state(raw: str) -> str:
    s = raw.strip().upper()
    # "CANCELLED by 12345" -> "CANCELLED"
    s = s.split()[0] if s else s
    return STATE_BUCKETS.get(s, "OTHER")


def parse_states(text: str) -> list[tuple[str, str]]:
    """Return (jobid, state) pairs from sacct/squeue tabular text.

    Accepts both whitespace-delimited and ``-P`` pipe-delimited output. Skips
    header lines and ``.batch`` / ``.extern`` sub-steps so array tasks are
    counted once.
    """
    rows: list[tuple[str, str]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("|") if "|" in line else line.split()
        if len(parts) < 2:
            continue
        jobid = parts[0]
        # Header rows
        if jobid.upper() in {"JOBID", "JOB_ID"}:
            continue
        # Sub-steps like 12345.batch / 12345.extern / 12345.0
        if "." in jobid:
            continue
        # Find the state field: sacct -> parts[1]; squeue default -> ST is col 5.
        state_raw = None
        for cand in (parts[1], parts[-1]):
            bucket = _normalize_state(cand)
            if bucket != "OTHER" or cand.strip().upper() in STATE_BUCKETS:
                state_raw = cand
                break
        if state_raw is None:
            state_raw = parts[1]
        rows.append((jobid, _normalize_state(state_raw)))
    return rows


def summarize(text: str) -> Counter:
    counts = Counter()
    for _, bucket in parse_states(text):
        counts[bucket] += 1
    return counts


def format_summary(counts: Counter) -> str:
    total = sum(counts.values())
    lines = [f"Array tasks: {total}"]
    for bucket in BUCKET_ORDER:
        if counts.get(bucket):
            pct = 100.0 * counts[bucket] / total if total else 0.0
            lines.append(f"  {bucket:10s} {counts[bucket]:5d}  ({pct:5.1f}%)")
    if total and counts.get("COMPLETED", 0) + counts.get("FAILED", 0) == total:
        lines.append("All tasks finished.")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Summarize SLURM array-job progress from sacct/squeue output."
    )
    ap.add_argument(
        "input", nargs="?", default="-",
        help="sacct/squeue text file, or '-' for stdin",
    )
    args = ap.parse_args(argv)

    text = sys.stdin.read() if args.input == "-" else open(args.input).read()
    print(format_summary(summarize(text)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
