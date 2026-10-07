from __future__ import annotations

from pathlib import Path
import csv
import json
import random
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bst_id.region_algebra import (  # noqa: E402
    AXES,
    Cell,
    compatible,
    contains,
    intersection,
    subsumes,
)
from bst_id.region_index import (  # noqa: E402
    RegionPrefixIndex,
    intersection_indexed,
)


SEED = 261007
RNG = random.Random(SEED)
RESULTS = []
FAILURES = []
PAIR_STATS = []


def record(name, cases, passed, seconds):
    RESULTS.append(
        {
            "test": name,
            "cases": cases,
            "passed": passed,
            "failed": cases - passed,
            "seconds": round(seconds, 6),
        }
    )


def rc(dims, zmin=1, zmax=10):
    parts = {}
    for axis in dims:
        z = RNG.randint(zmin, zmax)
        parts[axis] = (z, RNG.randrange(1 << z))
    return Cell.from_parts(**parts)


def random_region(dims, n, zmin=1, zmax=10):
    out = set()
    while len(out) < n:
        out.add(rc(dims, zmin, zmax))
    return tuple(out)


def fail(tag, *payload):
    if len(FAILURES) < 50:
        FAILURES.append((tag,) + payload)


# ---------------------------------------------------------------------------
# 1. Indexed ancestor / descendant / compatibility queries vs linear reference
# ---------------------------------------------------------------------------

for dims, regions, queries_per_region in [
    (("x",), 60, 50),
    (("x", "y"), 60, 50),
    (("x", "y", "f"), 60, 50),
    (("x", "y", "f", "t"), 60, 50),
]:
    t0 = time.perf_counter()
    cases = passed = 0

    for _ in range(regions):
        region = random_region(dims, RNG.randint(8, 40), 1, 10)
        index = RegionPrefixIndex(region)

        for _ in range(queries_per_region):
            q = rc(dims, 1, 12)

            got = set(index.subsuming_cells(q))
            exp = {r for r in region if subsumes(r, q)}
            cases += 1
            ok = got == exp
            passed += ok
            if not ok:
                fail(f"{len(dims)}D subsuming", q, got, exp)

            got = set(index.subsumed_cells(q))
            exp = {r for r in region if subsumes(q, r)}
            cases += 1
            ok = got == exp
            passed += ok
            if not ok:
                fail(f"{len(dims)}D subsumed", q, got, exp)

            got = set(index.compatible_candidates(q))
            exp = {r for r in region if compatible(r, q)}
            cases += 1
            ok = got == exp
            passed += ok
            if not ok:
                fail(f"{len(dims)}D compatible", q, got, exp)

            got = index.contains(q)
            exp = contains(region, q)
            cases += 1
            ok = got == exp
            passed += ok
            if not ok:
                fail(f"{len(dims)}D contains", q, got, exp)

    record(
        f"{len(dims)}D indexed prefix queries vs linear reference",
        cases,
        passed,
        time.perf_counter() - t0,
    )


# ---------------------------------------------------------------------------
# 2. Indexed Intersection vs current all-pairs reference implementation
# ---------------------------------------------------------------------------

for dims, trials in [
    (("x",), 600),
    (("x", "y"), 600),
    (("x", "y", "f"), 600),
    (("x", "y", "f", "t"), 600),
]:
    t0 = time.perf_counter()
    cases = passed = 0
    naive_pairs_total = 0
    indexed_pairs_total = 0

    for _ in range(trials):
        na = RNG.randint(2, 30)
        nb = RNG.randint(2, 30)
        a = random_region(dims, na, 1, 10)
        b = random_region(dims, nb, 1, 10)

        reference = intersection(a, b)
        indexed = intersection_indexed(a, b)

        cases += 1
        ok = indexed == reference
        passed += ok
        if not ok:
            fail(f"{len(dims)}D intersection", a, b, indexed, reference)

        # Candidate-pair count only; this is not a wall-clock speedup claim.
        if len(a) <= len(b):
            idx = RegionPrefixIndex(a)
            indexed_pairs = idx.count_compatible_pairs(b)
        else:
            idx = RegionPrefixIndex(b)
            indexed_pairs = idx.count_compatible_pairs(a)

        naive_pairs = len(a) * len(b)
        naive_pairs_total += naive_pairs
        indexed_pairs_total += indexed_pairs

    fraction = (
        indexed_pairs_total / naive_pairs_total if naive_pairs_total else 0.0
    )
    PAIR_STATS.append(
        {
            "dimensions": len(dims),
            "trials": trials,
            "naive_pairs": naive_pairs_total,
            "indexed_candidate_pairs": indexed_pairs_total,
            "candidate_fraction": fraction,
            "candidate_reduction_fraction": 1.0 - fraction,
        }
    )

    record(
        f"{len(dims)}D indexed Intersection vs all-pairs reference",
        cases,
        passed,
        time.perf_counter() - t0,
    )


# ---------------------------------------------------------------------------
# 3. Incremental add/discard integrity
# ---------------------------------------------------------------------------

for dims in [
    ("x",),
    ("x", "y"),
    ("x", "y", "f"),
    ("x", "y", "f", "t"),
]:
    t0 = time.perf_counter()
    cases = passed = 0
    pool = list(random_region(dims, 200, 1, 10))
    active = set(pool[:80])
    index = RegionPrefixIndex(active)

    for step in range(500):
        if step % 2 == 0:
            c = pool[RNG.randrange(len(pool))]
            active.add(c)
            index.add(c)
        else:
            if active:
                c = tuple(active)[RNG.randrange(len(active))]
                active.discard(c)
                index.discard(c)

        q = rc(dims, 1, 12)

        checks = [
            index.cells == frozenset(active),
            index.contains(q) == contains(active, q),
            set(index.subsuming_cells(q))
            == {r for r in active if subsumes(r, q)},
            set(index.compatible_candidates(q))
            == {r for r in active if compatible(r, q)},
        ]

        cases += len(checks)
        passed += sum(checks)
        if not all(checks):
            fail(f"{len(dims)}D update", step, q, checks)

    record(
        f"{len(dims)}D incremental index add/discard integrity",
        cases,
        passed,
        time.perf_counter() - t0,
    )


# ---------------------------------------------------------------------------
# Save outputs separately from the existing 38,732-case journal validation.
# ---------------------------------------------------------------------------

OUTDIR = Path(__file__).resolve().parent
summary = {
    "seed": SEED,
    "results": RESULTS,
    "pair_stats": PAIR_STATS,
    "failures": FAILURES,
}

with open(OUTDIR / "region_index_validation.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, default=str)

with open(
    OUTDIR / "region_index_validation_summary.csv",
    "w",
    newline="",
    encoding="utf-8",
) as f:
    writer = csv.DictWriter(
        f, fieldnames=["test", "cases", "passed", "failed", "seconds"]
    )
    writer.writeheader()
    writer.writerows(RESULTS)

with open(
    OUTDIR / "region_index_candidate_pairs.csv",
    "w",
    newline="",
    encoding="utf-8",
) as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "dimensions",
            "trials",
            "naive_pairs",
            "indexed_candidate_pairs",
            "candidate_fraction",
            "candidate_reduction_fraction",
        ],
    )
    writer.writeheader()
    writer.writerows(PAIR_STATS)

print(json.dumps(RESULTS, indent=2))
print()
print("Candidate-pair statistics:")
print(json.dumps(PAIR_STATS, indent=2))
print()
print(
    "TOTAL",
    sum(r["cases"] for r in RESULTS),
    "PASSED",
    sum(r["passed"] for r in RESULTS),
    "FAILED",
    sum(r["failed"] for r in RESULTS),
)
print("failures saved:", len(FAILURES))

if FAILURES:
    sys.exit(1)
