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
    Cell,
    compatible,
    contains,
    dense_difference_oracle,
    dense_intersection_oracle,
    difference,
    intersection,
    normalize,
    region_atoms,
    subsumes,
    union,
)


AXES = ("x", "y", "f", "t")
SEED = 261015
RESULTS = []
FAILURES = []


def record(name, cases, passed, seconds, oracle):
    RESULTS.append(
        {
            "test": name,
            "cases": cases,
            "passed": passed,
            "failed": cases - passed,
            "seconds": round(seconds, 6),
            "reference_oracle": oracle,
        }
    )


def fail(tag, *payload):
    if len(FAILURES) < 50:
        FAILURES.append((tag,) + payload)


def rc(rng, dims, zmin=2, zmax=4):
    kw = {}
    for axis in dims:
        z = rng.randint(zmin, zmax)
        kw[axis] = (z, rng.randrange(1 << z))
    return Cell.from_parts(**kw)


def random_region(rng, dims, n, zmin=2, zmax=4):
    return [rc(rng, dims, zmin, zmax) for _ in range(n)]


def union_oracle(a, b, w):
    return region_atoms(a, w) | region_atoms(b, w)


def anisotropic_pair(rng):
    x_parent = rng.randrange(1 << 2)
    y_parent = rng.randrange(1 << 2)
    x_suffix = rng.randrange(1 << 2)
    y_suffix = rng.randrange(1 << 2)

    a = Cell.from_parts(
        x=(4, (x_parent << 2) | x_suffix),
        y=(2, y_parent),
    )
    b = Cell.from_parts(
        x=(2, x_parent),
        y=(4, (y_parent << 2) | y_suffix),
    )
    return a, b


# ---------------------------------------------------------------------------
# R1-5a: multi-cell Boolean operands
# ---------------------------------------------------------------------------

rng = random.Random(SEED)
t0 = time.perf_counter()
cases = passed = 0
dims = ("x", "y")
w = {"x": 5, "y": 5}

for _ in range(600):
    a = random_region(rng, dims, rng.randint(2, 6))
    b = random_region(rng, dims, rng.randint(2, 6))

    checks = [
        (
            region_atoms(union(a, b), w),
            union_oracle(a, b, w),
            "union",
        ),
        (
            region_atoms(intersection(a, b), w),
            dense_intersection_oracle(a, b, w),
            "intersection",
        ),
        (
            region_atoms(difference(a, b, working_zoom=w), w),
            dense_difference_oracle(a, b, w),
            "difference",
        ),
    ]

    for got, expected, tag in checks:
        cases += 1
        ok = got == expected
        passed += ok
        if not ok:
            fail("multi-cell Boolean", tag, a, b, got, expected)

record(
    "2D multi-cell Boolean operands",
    cases,
    passed,
    time.perf_counter() - t0,
    "Dense Boolean atom sets",
)


# ---------------------------------------------------------------------------
# R1-5b: multiple subtractors
# ---------------------------------------------------------------------------

rng = random.Random(SEED + 1)
t0 = time.perf_counter()
cases = passed = 0
dims = ("x", "y", "f")
w = {"x": 4, "y": 4, "f": 4}

for _ in range(600):
    a = random_region(rng, dims, rng.randint(2, 6), 2, 4)
    b = random_region(rng, dims, rng.randint(2, 5), 2, 4)

    got = region_atoms(difference(a, b, working_zoom=w), w)
    expected = region_atoms(a, w) - region_atoms(b, w)

    cases += 1
    ok = got == expected
    passed += ok
    if not ok:
        fail("multiple subtractors", a, b, got, expected)

record(
    "3D Difference with multiple subtractors",
    cases,
    passed,
    time.perf_counter() - t0,
    "Dense set difference",
)


# ---------------------------------------------------------------------------
# R1-5c: cross-nested anisotropic overlap with multi-cell operands
# ---------------------------------------------------------------------------

rng = random.Random(SEED + 2)
t0 = time.perf_counter()
cases = passed = 0
w = {"x": 5, "y": 5}

for _ in range(600):
    a1, b1 = anisotropic_pair(rng)
    a2, b2 = anisotropic_pair(rng)
    a = [a1, a2]
    b = [b1, b2]

    construction_ok = (
        compatible(a1, b1)
        and not subsumes(a1, b1)
        and not subsumes(b1, a1)
    )

    checks = [
        (
            construction_ok
            and region_atoms(intersection(a, b), w)
            == dense_intersection_oracle(a, b, w),
            "intersection",
        ),
        (
            construction_ok
            and region_atoms(difference(a, b, working_zoom=w), w)
            == dense_difference_oracle(a, b, w),
            "difference",
        ),
    ]

    for ok, tag in checks:
        cases += 1
        passed += ok
        if not ok:
            fail("anisotropic overlap", tag, a, b)

record(
    "2D anisotropic cross-nested multi-cell overlap",
    cases,
    passed,
    time.perf_counter() - t0,
    "Dense Boolean atom sets + cross-nesting checks",
)


# ---------------------------------------------------------------------------
# R1-5d: chained region operations
# ---------------------------------------------------------------------------

rng = random.Random(SEED + 3)
t0 = time.perf_counter()
cases = passed = 0
dims = ("x", "y", "f")
w = {"x": 4, "y": 4, "f": 4}

for _ in range(600):
    a = random_region(rng, dims, rng.randint(2, 5), 2, 4)
    b = random_region(rng, dims, rng.randint(2, 5), 2, 4)
    c = random_region(rng, dims, rng.randint(2, 4), 2, 4)
    d = random_region(rng, dims, rng.randint(2, 5), 2, 4)

    stage1 = union(a, b)
    stage2 = difference(stage1, c, working_zoom=w)
    stage3 = intersection(stage2, d)
    final = normalize(stage3)

    got = region_atoms(final, w)
    expected = (
        (region_atoms(a, w) | region_atoms(b, w))
        - region_atoms(c, w)
    ) & region_atoms(d, w)

    ok = got == expected and normalize(final) == final
    cases += 1
    passed += ok
    if not ok:
        fail("chained operations", a, b, c, d, final, expected)

record(
    "3D chained Union-Difference-Intersection-Normalize",
    cases,
    passed,
    time.perf_counter() - t0,
    "Composed dense atom-set semantics",
)


# ---------------------------------------------------------------------------
# R1-5e: UAS negative conformance cases
# ---------------------------------------------------------------------------

reserved = [
    Cell.from_parts(
        x=(3, 0b101),
        y=(4, 0b0110),
        f=(3, 0b010),
        t=(3, 0b011),
    ),
    Cell.from_parts(
        x=(3, 0b001),
        y=(4, 0b1100),
        f=(3, 0b110),
        t=(3, 0b110),
    ),
]

# XY and T inside, height outside.
t0 = time.perf_counter()
cases = passed = 0
for sx in range(4):
    for sy in range(4):
        for sf in range(4):
            for st in range(4):
                q = Cell.from_parts(
                    x=(5, (0b101 << 2) | sx),
                    y=(6, (0b0110 << 2) | sy),
                    f=(5, (0b011 << 2) | sf),
                    t=(5, (0b011 << 2) | st),
                )
                cases += 1
                ok = not contains(reserved, q)
                passed += ok
                if not ok:
                    fail("XY-in/H-out", q)

record(
    "XYHT conformance negative: XY and T inside, H outside",
    cases,
    passed,
    time.perf_counter() - t0,
    "Binary subsumption (expected False)",
)

# XYH inside, time outside.
t0 = time.perf_counter()
cases = passed = 0
for sx in range(4):
    for sy in range(4):
        for sf in range(4):
            for st in range(4):
                q = Cell.from_parts(
                    x=(5, (0b101 << 2) | sx),
                    y=(6, (0b0110 << 2) | sy),
                    f=(5, (0b010 << 2) | sf),
                    t=(5, (0b100 << 2) | st),
                )
                cases += 1
                ok = not contains(reserved, q)
                passed += ok
                if not ok:
                    fail("XYH-in/T-out", q)

record(
    "XYHT conformance negative: XYH inside, T outside",
    cases,
    passed,
    time.perf_counter() - t0,
    "Binary subsumption (expected False)",
)


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

OUTDIR = Path(__file__).resolve().parent
payload = {
    "seed": SEED,
    "reviewer": "Reviewer 1, Comment 5",
    "results": RESULTS,
    "failures": FAILURES,
}

with open(
    OUTDIR / "reviewer1_5_validation.json",
    "w",
    encoding="utf-8",
) as f:
    json.dump(payload, f, indent=2, default=str)

with open(
    OUTDIR / "reviewer1_5_validation_summary.csv",
    "w",
    newline="",
    encoding="utf-8",
) as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "test",
            "cases",
            "passed",
            "failed",
            "seconds",
            "reference_oracle",
        ],
    )
    writer.writeheader()
    writer.writerows(RESULTS)

print(json.dumps(RESULTS, indent=2))
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
