from __future__ import annotations

import random

from bst_id.region_algebra import (
    Cell,
    compatible,
    contains,
    intersection,
    subsumes,
)
from bst_id.region_index import RegionPrefixIndex, intersection_indexed


AXES = ("x", "y", "f", "t")


def rc(rng: random.Random, dims, zmin=1, zmax=8) -> Cell:
    parts = {}
    for axis in dims:
        z = rng.randint(zmin, zmax)
        parts[axis] = (z, rng.randrange(1 << z))
    return Cell.from_parts(**parts)


def test_indexed_subsumption_matches_linear_scan():
    region = [
        Cell.from_parts(x=(2, 0b10), y=(2, 0b01)),
        Cell.from_parts(x=(3, 0b101), y=(3, 0b011)),
        Cell.from_parts(x=(4, 0b1110), y=(2, 0b10)),
    ]
    index = RegionPrefixIndex(region)

    queries = [
        Cell.from_parts(x=(5, 0b10111), y=(5, 0b01101)),
        Cell.from_parts(x=(4, 0b1110), y=(4, 0b1011)),
        Cell.from_parts(x=(4, 0b0001), y=(4, 0b0001)),
    ]

    for q in queries:
        expected = {r for r in region if subsumes(r, q)}
        assert set(index.subsuming_cells(q)) == expected
        assert index.contains(q) == contains(region, q)


def test_indexed_descendants_match_linear_scan():
    region = [
        Cell.from_parts(x=(2, 0b10), y=(2, 0b01)),
        Cell.from_parts(x=(3, 0b101), y=(3, 0b011)),
        Cell.from_parts(x=(4, 0b1011), y=(4, 0b0110)),
        Cell.from_parts(x=(3, 0b111), y=(3, 0b010)),
    ]
    index = RegionPrefixIndex(region)
    container = Cell.from_parts(x=(2, 0b10), y=(2, 0b01))

    expected = {r for r in region if subsumes(container, r)}
    assert set(index.subsumed_cells(container)) == expected


def test_compatible_candidates_match_linear_scan_anisotropic():
    region = [
        Cell.from_parts(x=(3, 0b101), y=(2, 0b10)),
        Cell.from_parts(x=(2, 0b10), y=(4, 0b1011)),
        Cell.from_parts(x=(3, 0b111), y=(3, 0b000)),
        Cell.from_parts(x=(4, 0b1011), y=(4, 0b1010)),
    ]
    index = RegionPrefixIndex(region)
    query = Cell.from_parts(x=(4, 0b1011), y=(4, 0b1011))

    expected = {r for r in region if compatible(r, query)}
    assert set(index.compatible_candidates(query)) == expected


def test_indexed_intersection_matches_reference():
    a = [
        Cell.from_parts(x=(3, 0b101), y=(2, 0b10)),
        Cell.from_parts(x=(3, 0b111), y=(3, 0b000)),
    ]
    b = [
        Cell.from_parts(x=(2, 0b10), y=(4, 0b1011)),
        Cell.from_parts(x=(4, 0b1110), y=(2, 0b00)),
    ]

    assert intersection_indexed(a, b) == intersection(a, b)


def test_add_discard_update():
    a = Cell.from_parts(x=(2, 0b10), y=(2, 0b01))
    b = Cell.from_parts(x=(3, 0b101), y=(3, 0b011))
    q = Cell.from_parts(x=(5, 0b10111), y=(5, 0b01101))

    index = RegionPrefixIndex([a])
    assert index.contains(q)

    index.add(b)
    assert len(index) == 2

    index.discard(a)
    assert len(index) == 1
    assert index.contains(q)

    index.discard(b)
    assert len(index) == 0
    assert not index.contains(q)


def test_randomized_index_queries_match_reference():
    rng = random.Random(261007)

    for dims in [
        ("x",),
        ("x", "y"),
        ("x", "y", "f"),
        ("x", "y", "f", "t"),
    ]:
        for _ in range(40):
            region = {rc(rng, dims, 1, 7) for _ in range(24)}
            index = RegionPrefixIndex(region)

            for _ in range(20):
                q = rc(rng, dims, 1, 8)

                assert set(index.subsuming_cells(q)) == {
                    r for r in region if subsumes(r, q)
                }
                assert set(index.subsumed_cells(q)) == {
                    r for r in region if subsumes(q, r)
                }
                assert set(index.compatible_candidates(q)) == {
                    r for r in region if compatible(r, q)
                }
                assert index.contains(q) == contains(region, q)


def test_randomized_indexed_intersection_matches_reference():
    rng = random.Random(261008)

    for dims in [
        ("x",),
        ("x", "y"),
        ("x", "y", "f"),
        ("x", "y", "f", "t"),
    ]:
        for _ in range(100):
            a = {rc(rng, dims, 1, 7) for _ in range(rng.randint(1, 12))}
            b = {rc(rng, dims, 1, 7) for _ in range(rng.randint(1, 12))}
            assert intersection_indexed(a, b) == intersection(a, b)


def test_presence_domain_mismatch_rejected():
    index = RegionPrefixIndex([Cell.from_parts(x=(3, 0b101))])
    bad = Cell.from_parts(x=(3, 0b101), y=(3, 0b010))

    try:
        index.contains(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("presence-domain mismatch must raise ValueError")
