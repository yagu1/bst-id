from __future__ import annotations

import random

from bst_id.region_algebra import Cell, normalize
from bst_id.region_index import (
    normalize_indexed,
    normalize_indexed_batch,
)


def rc(rng, dims, zmin=1, zmax=7):
    kw = {}
    for axis in dims:
        z = rng.randint(zmin, zmax)
        kw[axis] = (z, rng.randrange(1 << z))
    return Cell.from_parts(**kw)


def test_empty_and_singleton():
    assert normalize_indexed([]) == normalize([])
    assert normalize_indexed_batch([]) == normalize([])

    c = Cell.from_parts(x=(3, 0b101))
    assert normalize_indexed([c]) == normalize([c])
    assert normalize_indexed_batch([c]) == normalize([c])


def test_subsumption_case_matches_reference():
    cells = [
        Cell.from_parts(x=(2, 0b10)),
        Cell.from_parts(x=(3, 0b100)),
        Cell.from_parts(x=(4, 0b1010)),
        Cell.from_parts(x=(5, 0b10111)),
    ]
    ref = normalize(cells)
    assert normalize_indexed(cells) == ref
    assert normalize_indexed_batch(cells) == ref


def test_merge_chain_matches_reference():
    cells = [Cell.from_parts(x=(5, i)) for i in range(32)]
    ref = normalize(cells)
    assert normalize_indexed(cells) == ref
    assert normalize_indexed_batch(cells) == ref


def test_anisotropic_merge_matches_reference():
    cells = [
        Cell.from_parts(x=(4, x), y=(3, 0b101))
        for x in range(16)
    ]
    ref = normalize(cells)
    assert normalize_indexed(cells) == ref
    assert normalize_indexed_batch(cells) == ref


def test_randomized_exact_output_matches_reference():
    rng = random.Random(261009)

    for dims in [
        ("x",),
        ("x", "y"),
        ("x", "y", "f"),
        ("x", "y", "f", "t"),
    ]:
        for _ in range(500):
            cells = [
                rc(rng, dims, 1, 7)
                for _ in range(rng.randint(1, 30))
            ]

            ref = normalize(cells)
            got = normalize_indexed(cells)
            batch = normalize_indexed_batch(cells)

            assert got == ref
            assert batch == ref


def test_input_order_invariance():
    rng = random.Random(261010)
    dims = ("x", "y", "f", "t")
    cells = [rc(rng, dims, 1, 6) for _ in range(24)]

    ref = normalize(cells)
    assert normalize_indexed(cells) == ref
    assert normalize_indexed(list(reversed(cells))) == ref
    assert normalize_indexed_batch(cells) == ref
    assert normalize_indexed_batch(list(reversed(cells))) == ref
