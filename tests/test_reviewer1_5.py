from __future__ import annotations

import random

from bst_id.region_algebra import (
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


def _rc(rng, dims, zmin=2, zmax=4):
    parts = {}
    for axis in dims:
        z = rng.randint(zmin, zmax)
        parts[axis] = (z, rng.randrange(1 << z))
    return Cell.from_parts(**parts)


def _random_region(rng, dims, n, zmin=2, zmax=4):
    return [_rc(rng, dims, zmin, zmax) for _ in range(n)]


def _union_oracle(a, b, working_zoom):
    return region_atoms(a, working_zoom) | region_atoms(b, working_zoom)


def test_multicell_boolean_operands_against_dense_oracle():
    rng = random.Random(261015)
    dims = ("x", "y")
    wz = {"x": 5, "y": 5}

    for _ in range(100):
        a = _random_region(rng, dims, rng.randint(2, 6))
        b = _random_region(rng, dims, rng.randint(2, 6))

        assert region_atoms(union(a, b), wz) == _union_oracle(a, b, wz)
        assert region_atoms(intersection(a, b), wz) == dense_intersection_oracle(a, b, wz)
        assert region_atoms(difference(a, b, working_zoom=wz), wz) == dense_difference_oracle(a, b, wz)


def test_multiple_subtractors_against_dense_oracle():
    rng = random.Random(261016)
    dims = ("x", "y", "f")
    wz = {"x": 4, "y": 4, "f": 4}

    for _ in range(100):
        minuend = _random_region(rng, dims, rng.randint(2, 6), 2, 4)
        subtractors = _random_region(rng, dims, rng.randint(2, 5), 2, 4)

        got = region_atoms(
            difference(minuend, subtractors, working_zoom=wz),
            wz,
        )
        expected = region_atoms(minuend, wz) - region_atoms(subtractors, wz)
        assert got == expected


def _anisotropic_pair(rng):
    # a is finer in x and coarser in y; b is the opposite.
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


def test_anisotropic_cross_nested_multicell_overlap():
    rng = random.Random(261017)
    wz = {"x": 5, "y": 5}

    for _ in range(100):
        a1, b1 = _anisotropic_pair(rng)
        a2, b2 = _anisotropic_pair(rng)
        a = [a1, a2]
        b = [b1, b2]

        # At least the constructed paired cells are compatible and cross-nested.
        assert compatible(a1, b1)
        assert not subsumes(a1, b1)
        assert not subsumes(b1, a1)

        assert region_atoms(intersection(a, b), wz) == dense_intersection_oracle(a, b, wz)
        assert region_atoms(difference(a, b, working_zoom=wz), wz) == dense_difference_oracle(a, b, wz)


def test_chained_operations_against_dense_oracle():
    rng = random.Random(261018)
    dims = ("x", "y", "f")
    wz = {"x": 4, "y": 4, "f": 4}

    for _ in range(100):
        a = _random_region(rng, dims, rng.randint(2, 5), 2, 4)
        b = _random_region(rng, dims, rng.randint(2, 5), 2, 4)
        c = _random_region(rng, dims, rng.randint(2, 4), 2, 4)
        d = _random_region(rng, dims, rng.randint(2, 5), 2, 4)

        # Union -> Difference -> Intersection -> Normalize
        stage1 = union(a, b)
        stage2 = difference(stage1, c, working_zoom=wz)
        stage3 = intersection(stage2, d)
        got = region_atoms(normalize(stage3), wz)

        expected = (
            (region_atoms(a, wz) | region_atoms(b, wz))
            - region_atoms(c, wz)
        ) & region_atoms(d, wz)

        assert got == expected
        assert normalize(stage3) == stage3


def test_uas_xy_inside_but_height_outside():
    reserved = [
        Cell.from_parts(
            x=(3, 0b101),
            y=(4, 0b0110),
            f=(3, 0b010),
            t=(3, 0b011),
        ),
        # A second far-away reservation verifies multi-cell membership does not
        # accidentally accept the negative query.
        Cell.from_parts(
            x=(3, 0b001),
            y=(4, 0b1100),
            f=(3, 0b110),
            t=(3, 0b110),
        ),
    ]

    for sx in range(4):
        for sy in range(4):
            for sf in range(4):
                for st in range(4):
                    q = Cell.from_parts(
                        x=(5, (0b101 << 2) | sx),       # inside X
                        y=(6, (0b0110 << 2) | sy),      # inside Y
                        f=(5, (0b011 << 2) | sf),       # outside H/F
                        t=(5, (0b011 << 2) | st),       # inside T
                    )
                    assert not contains(reserved, q)


def test_uas_xyh_inside_but_time_outside():
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

    for sx in range(4):
        for sy in range(4):
            for sf in range(4):
                for st in range(4):
                    q = Cell.from_parts(
                        x=(5, (0b101 << 2) | sx),       # inside X
                        y=(6, (0b0110 << 2) | sy),      # inside Y
                        f=(5, (0b010 << 2) | sf),       # inside H/F
                        t=(5, (0b100 << 2) | st),       # outside T
                    )
                    assert not contains(reserved, q)
