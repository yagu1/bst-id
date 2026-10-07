from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bst_id.region_algebra import Cell, normalize  # noqa: E402
from bst_id.region_index import (  # noqa: E402
    normalize_indexed_profile,
    normalize_indexed_batch_profile,
)


def antichain_1d(n):
    # Same depth, even values only: no duplicates, no subsumption, no sibling pairs.
    z = max(2, math.ceil(math.log2(max(2, 2 * n))))
    return [Cell.from_parts(x=(z, 2 * i)) for i in range(n)]


def antichain_4d(n):
    # Vary x only and keep the other axes fixed. Even x values avoid sibling pairs.
    z = max(2, math.ceil(math.log2(max(2, 2 * n))))
    return [
        Cell.from_parts(
            x=(z, 2 * i),
            y=(8, 0b01010101),
            f=(8, 0b00110011),
            t=(8, 0b00011100),
        )
        for i in range(n)
    ]


def merge_heavy_1d(n):
    # Use the largest power-of-two complete sibling block <= n.
    p = 1 << int(math.floor(math.log2(max(2, n))))
    z = int(math.log2(p))
    z = max(1, z)
    return [Cell.from_parts(x=(z, i)) for i in range(1 << z)]


def anisotropic_merge_2d(n):
    p = 1 << int(math.floor(math.log2(max(2, n))))
    z = int(math.log2(p))
    z = max(1, z)
    return [
        Cell.from_parts(x=(z, i), y=(7, 0b0101010))
        for i in range(1 << z)
    ]


def mixed_random_4d(n, seed):
    rng = random.Random(seed)
    out = set()
    while len(out) < n:
        kw = {}
        for axis in ("x", "y", "f", "t"):
            z = rng.randint(3, 10)
            kw[axis] = (z, rng.randrange(1 << z))
        out.add(Cell.from_parts(**kw))
    return list(out)


SCENARIOS = {
    "antichain_1d": antichain_1d,
    "antichain_4d": antichain_4d,
    "merge_heavy_1d": merge_heavy_1d,
    "anisotropic_merge_2d": anisotropic_merge_2d,
}


def timed(fn, cells):
    t0 = time.perf_counter()
    out = fn(cells)
    return out, time.perf_counter() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", nargs="+", type=int, default=[32, 64, 128, 256])
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--include-random", action="store_true")
    args = ap.parse_args()

    outdir = Path(__file__).resolve().parent
    raw_rows = []
    summaries = []

    scenario_items = list(SCENARIOS.items())
    if args.include_random:
        scenario_items.append(
            ("mixed_random_4d", lambda n: mixed_random_4d(n, 261011 + n))
        )

    for scenario, make_cells in scenario_items:
        for requested_n in args.sizes:
            cells = make_cells(requested_n)
            n = len(set(cells))

            # Compute reference once outside the timed repeats for correctness anchor.
            reference = normalize(cells)

            algorithms = [
                ("reference", lambda xs: (normalize(xs), None)),
                ("indexed_n1", normalize_indexed_profile),
                ("indexed_n1_batch_n2", normalize_indexed_batch_profile),
            ]

            for name, fn in algorithms:
                times = []
                exact_flags = []
                stats_last = None

                for rep in range(args.repeats):
                    t0 = time.perf_counter()
                    out, stats = fn(cells)
                    dt = time.perf_counter() - t0

                    times.append(dt)
                    exact = out == reference
                    exact_flags.append(exact)
                    stats_last = stats

                    raw_rows.append(
                        {
                            "scenario": scenario,
                            "requested_n": requested_n,
                            "input_n": n,
                            "algorithm": name,
                            "repeat": rep,
                            "seconds": dt,
                            "output_n": len(out),
                            "exact_equal_reference": exact,
                            "passes": "" if stats is None else stats.passes,
                            "subsumption_queries": ""
                            if stats is None
                            else stats.subsumption_queries,
                            "subsumed_removed": ""
                            if stats is None
                            else stats.subsumed_removed,
                            "index_insertions": ""
                            if stats is None
                            else stats.index_insertions,
                            "sibling_merges": ""
                            if stats is None
                            else stats.sibling_merges,
                        }
                    )

                summaries.append(
                    {
                        "scenario": scenario,
                        "requested_n": requested_n,
                        "input_n": n,
                        "algorithm": name,
                        "median_seconds": statistics.median(times),
                        "min_seconds": min(times),
                        "max_seconds": max(times),
                        "output_n": len(reference),
                        "all_exact_equal_reference": all(exact_flags),
                        "passes": ""
                        if stats_last is None
                        else stats_last.passes,
                        "subsumption_queries": ""
                        if stats_last is None
                        else stats_last.subsumption_queries,
                        "subsumed_removed": ""
                        if stats_last is None
                        else stats_last.subsumed_removed,
                        "index_insertions": ""
                        if stats_last is None
                        else stats_last.index_insertions,
                        "sibling_merges": ""
                        if stats_last is None
                        else stats_last.sibling_merges,
                    }
                )

    # Add speed ratio relative to reference median for each scenario/size.
    ref_time = {
        (r["scenario"], r["requested_n"]): r["median_seconds"]
        for r in summaries
        if r["algorithm"] == "reference"
    }
    for row in summaries:
        base = ref_time[(row["scenario"], row["requested_n"])]
        row["reference_time_over_algorithm_time"] = (
            base / row["median_seconds"] if row["median_seconds"] > 0 else ""
        )

    raw_path = outdir / "normalize_benchmark_raw.csv"
    summary_path = outdir / "normalize_benchmark_summary.csv"
    json_path = outdir / "normalize_benchmark.json"

    with open(raw_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=raw_rows[0].keys())
        writer.writeheader()
        writer.writerows(raw_rows)

    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=summaries[0].keys())
        writer.writeheader()
        writer.writerows(summaries)

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(
            {
                "sizes": args.sizes,
                "repeats": args.repeats,
                "summary": summaries,
            },
            f,
            indent=2,
        )

    print("Wrote:")
    print(" ", raw_path)
    print(" ", summary_path)
    print(" ", json_path)
    print()

    failed = [r for r in summaries if not r["all_exact_equal_reference"]]
    if failed:
        print("ERROR: indexed normalization differed from reference:")
        print(json.dumps(failed, indent=2))
        raise SystemExit(1)

    # Compact console summary.
    for row in summaries:
        print(
            f'{row["scenario"]:24s} '
            f'N={row["input_n"]:4d} '
            f'{row["algorithm"]:20s} '
            f'{row["median_seconds"]:.6f}s '
            f'x{row["reference_time_over_algorithm_time"]:.3f}'
        )


if __name__ == "__main__":
    main()
