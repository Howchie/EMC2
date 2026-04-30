#!/usr/bin/env python3
"""Summarize recovery tables from EMC2 demo logs.

This script parses the `recovery(...)` output printed by the WorkingTests
demo scripts and reports, per test block:

* parameter truth / posterior median / 95% interval
* whether the true value fell outside the 95% interval
* block-level RMSE and coverage

It is intentionally lightweight so it can be reused on future demo logs.

Example:
    python tools/recovery_log_summary.py --log demo_cens_trunc.log \
        --script WorkingTests/demo_cens_trunc.R
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional


@dataclass
class RecoveryRow:
    block: int
    label: str
    param: str
    q25: float
    median: float
    q975: float
    true: float
    miss: float
    rmse: float
    coverage: float

    @property
    def outside95(self) -> bool:
        return not (self.q25 <= self.true <= self.q975)

    @property
    def abs_miss(self) -> float:
        return abs(self.miss)


def parse_labels(script_text: str) -> List[str]:
    """Extract quoted labels from `label = "..."` call sites in the demo script."""
    labels: List[str] = []
    for line in script_text.splitlines():
        s = line.strip()
        if not s.startswith("label = "):
            continue
        if "label = NULL" in s:
            continue
        if '"' not in s:
            continue
        value = s.split("=", 1)[1].strip().rstrip(",").strip()
        labels.append(value.strip('"'))
    return labels


def parse_recovery_blocks(log_text: str) -> List[tuple[list[str], dict[str, float]]]:
    """Return raw recovery blocks as (rows, stats) tuples."""
    lines = log_text.splitlines()
    blocks: List[tuple[list[str], dict[str, float]]] = []
    i = 0
    while i < len(lines):
        if lines[i].strip() != "$`1`$quantiles":
            i += 1
            continue

        j = i + 1
        rows: list[str] = []
        while j < len(lines) and lines[j].strip() != "$`1`$stats":
            if lines[j].strip():
                rows.append(lines[j])
            j += 1

        stats: dict[str, float] = {}
        if j + 2 < len(lines):
            stat_names = lines[j + 1].split()
            stat_vals = lines[j + 2].split()
            for name, value in zip(stat_names, stat_vals):
                try:
                    stats[name] = float(value)
                except ValueError:
                    continue

        blocks.append((rows, stats))
        i = j

    return blocks


def parse_rows(
    blocks: Iterable[tuple[list[str], dict[str, float]]],
    labels: Optional[List[str]] = None,
) -> List[RecoveryRow]:
    rows: List[RecoveryRow] = []
    for block_idx, (raw_rows, stats) in enumerate(blocks, 1):
        label = labels[block_idx - 1] if labels and block_idx - 1 < len(labels) else f"block_{block_idx:02d}"
        for raw in raw_rows[1:]:  # skip column header
            parts = raw.split()
            if len(parts) < 6:
                continue
            q25, median, q975, true, miss = map(float, parts[1:6])
            rows.append(
                RecoveryRow(
                    block=block_idx,
                    label=label,
                    param=parts[0],
                    q25=q25,
                    median=median,
                    q975=q975,
                    true=true,
                    miss=miss,
                    rmse=stats.get("rmse", float("nan")),
                    coverage=stats.get("coverage", float("nan")),
                )
            )
    return rows


def print_summary(rows: List[RecoveryRow]) -> None:
    by_block: dict[int, list[RecoveryRow]] = {}
    for row in rows:
        by_block.setdefault(row.block, []).append(row)

    print("block label outside_n max_abs_miss rmse coverage bad_params")
    for block in sorted(by_block):
        block_rows = by_block[block]
        label = block_rows[0].label
        outside = [r for r in block_rows if r.outside95]
        bad_params = ",".join(r.param for r in outside)
        max_abs_miss = max(r.abs_miss for r in block_rows)
        print(
            f"{block:>5} {label:<24} {len(outside):>9} "
            f"{max_abs_miss:>12.6f} {block_rows[0].rmse:>8.6f} "
            f"{block_rows[0].coverage:>8.3f} {bad_params}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", required=True, type=Path, help="Path to the demo log file")
    parser.add_argument(
        "--script",
        type=Path,
        help="Optional demo R script to source labels from, so blocks get meaningful names",
    )
    args = parser.parse_args()

    log_text = args.log.read_text()
    labels = parse_labels(args.script.read_text()) if args.script else None
    blocks = parse_recovery_blocks(log_text)
    rows = parse_rows(blocks, labels=labels)
    print_summary(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
