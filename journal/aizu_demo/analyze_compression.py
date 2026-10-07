# journal/aizu_demo/analyze_compression.py
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence, Set, Tuple

T_BITS = 32
DEFAULT_X_ZOOM = 21
DEFAULT_Y_ZOOM = 21
DEFAULT_H_ZOOM = 12

# A full XYHT prefix is represented as:
# ((zx, ix), (zy, iy), (zh, ih), (zt, it))
Prefix = Tuple[Tuple[int, int], Tuple[int, int], Tuple[int, int], Tuple[int, int]]
Attribute = Tuple[str, ...]
IntervalEntry = Tuple[int, int, str]


def bst_bits_3d(zx: int, zy: int, zh: int) -> int:
    """BST-ID payload/header bits for XYH (4 flags + 3*5 zoom bits + prefixes)."""
    return 4 + 3 * 5 + zx + zy + zh


def bst_bits_4d(zx: int, zy: int, zh: int, zt: int) -> int:
    """BST-ID payload/header bits for XYHT (4 flags + 4*5 zoom bits + prefixes)."""
    return 4 + 4 * 5 + zx + zy + zh + zt


def human_bits(bits: int) -> str:
    return f"{bits / 8.0 / 1024.0:.2f} KiB"


def dyadic_time_cover(start_s: int, end_s: int) -> List[Tuple[int, int]]:
    """Exact minimal-aligned dyadic cover of [start_s, end_s).

    Returns (zoom_t, t_index) pairs for the BST-ID 32-bit Unix-time hierarchy.
    A prefix at zoom z covers 2**(32-z) seconds.
    """
    if not (0 <= start_s < end_s <= (1 << T_BITS)):
        raise ValueError(f"time interval out of 32-bit Unix range: {start_s}, {end_s}")

    out: List[Tuple[int, int]] = []
    t = start_s
    while t < end_s:
        remaining = end_s - t
        largest_by_length = 1 << (remaining.bit_length() - 1)
        alignment = (t & -t) if t else (1 << T_BITS)
        block = min(largest_by_length, alignment)

        # The scenario intervals are short, so zoom 0 will never normally occur.
        # Keep the representation compatible with the BST-ID implementation
        # (zoom in [1,32]) by splitting the universe-size block if necessary.
        if block >= (1 << T_BITS):
            block = 1 << (T_BITS - 1)

        log2_block = block.bit_length() - 1
        zoom_t = T_BITS - log2_block
        if not (1 <= zoom_t <= 32):
            raise ValueError(f"invalid T zoom {zoom_t} for block {block}")
        t_index = t >> log2_block
        out.append((zoom_t, t_index))
        t += block

    return out


def prefix_duration_s(prefix: Prefix) -> int:
    zt = prefix[3][0]
    return 1 << (T_BITS - zt)


def sweep_attribute_segments(entries: Sequence[IntervalEntry]) -> Iterator[Tuple[int, int, Attribute]]:
    """Partition one fixed XYH atom into maximal intervals of equal flight-set attribute."""
    events: Dict[int, List[Tuple[int, str]]] = defaultdict(list)
    for start, end, flight_id in entries:
        if end <= start:
            continue
        events[start].append((+1, flight_id))
        events[end].append((-1, flight_id))

    active_counts: Counter[str] = Counter()
    previous_t: int | None = None
    previous_attr: Attribute | None = None
    previous_start: int | None = None

    def emit_pending(stop_t: int):
        nonlocal previous_attr, previous_start
        if previous_attr and previous_start is not None and stop_t > previous_start:
            result = (previous_start, stop_t, previous_attr)
            previous_attr = None
            previous_start = None
            return result
        return None

    for t in sorted(events):
        if previous_t is not None and t > previous_t:
            current_attr = tuple(sorted(fid for fid, count in active_counts.items() if count > 0))
            if current_attr:
                if previous_attr == current_attr:
                    # Extend the pending maximal segment.
                    pass
                else:
                    pending = emit_pending(previous_t)
                    if pending is not None:
                        yield pending
                    previous_attr = current_attr
                    previous_start = previous_t
            else:
                pending = emit_pending(previous_t)
                if pending is not None:
                    yield pending

        # Half-open interval semantics: update state at t for [t, next_t).
        for delta, flight_id in events[t]:
            active_counts[flight_id] += delta
            if active_counts[flight_id] <= 0:
                del active_counts[flight_id]

        previous_t = t

    if previous_t is not None:
        pending = emit_pending(previous_t)
        if pending is not None:
            yield pending


def merge_siblings_on_axis(cells: Set[Prefix], axis: int) -> Tuple[Set[Prefix], int]:
    """One binary sibling-merge layer on one axis.

    All mergeable sibling pairs for this axis are merged simultaneously.
    Inputs are a disjoint partition for one semantic attribute, so subsumption
    removal is unnecessary: sibling merging preserves a disjoint partition.
    """
    buckets: Dict[Prefix, Dict[int, Prefix]] = {}

    for cell in cells:
        zoom, value = cell[axis]
        if zoom <= 1:
            continue
        parent_desc = (zoom - 1, value >> 1)
        parent = list(cell)
        parent[axis] = parent_desc
        parent_key: Prefix = tuple(parent)  # type: ignore[assignment]
        buckets.setdefault(parent_key, {})[value & 1] = cell

    pairs: List[Tuple[Prefix, Prefix, Prefix]] = []
    for parent, bits in buckets.items():
        if 0 in bits and 1 in bits:
            pairs.append((bits[0], bits[1], parent))

    if not pairs:
        return cells, 0

    new_cells = set(cells)
    for c0, c1, parent in pairs:
        new_cells.discard(c0)
        new_cells.discard(c1)
        new_cells.add(parent)

    return new_cells, len(pairs)


def normalize_prefix_partition(cells: Set[Prefix]) -> Tuple[Set[Prefix], Dict[str, int]]:
    """Deterministic attribute-aware XYHT sibling normalization.

    Axis priority is x -> y -> h -> t.  After any merge layer, the algorithm
    restarts from x.  This is a deterministic normal form under the specified
    policy; it does not claim a globally minimum-cardinality rectangle cover.
    """
    axis_names = ("x", "y", "h", "t")
    merges = {name: 0 for name in axis_names}
    current = set(cells)

    while True:
        changed = False
        for axis, name in enumerate(axis_names):
            updated, count = merge_siblings_on_axis(current, axis)
            if count:
                current = updated
                merges[name] += count
                changed = True
                break
        if not changed:
            break

    return current, merges


def get_dynamic_entries(data: Mapping[str, object]) -> Tuple[
    Dict[Tuple[int, int, int], List[IntervalEntry]], int, int, int
]:
    """Collect per-flight dynamic atom intervals into a global fixed-grid index."""
    metadata = data.get("metadata", {})
    if not isinstance(metadata, Mapping):
        metadata = {}
    reservation_meta = metadata.get("reservation", {})
    if not isinstance(reservation_meta, Mapping):
        reservation_meta = {}

    x_zoom = int(reservation_meta.get("margin_xy_zoom", DEFAULT_X_ZOOM))
    y_zoom = x_zoom
    h_zoom = int(reservation_meta.get("h_zoom", DEFAULT_H_ZOOM))

    entries_by_atom: Dict[Tuple[int, int, int], List[IntervalEntry]] = defaultdict(list)

    flights = data.get("flights", [])
    if not isinstance(flights, list):
        raise ValueError("utm_scenario.json does not contain a flights list")

    for flight in flights:
        if not isinstance(flight, Mapping):
            continue
        flight_id = str(flight["flight_id"])
        reservation = flight.get("reservation", {})
        if not isinstance(reservation, Mapping):
            continue
        atoms = reservation.get("dynamic_atoms")
        if atoms is None:
            raise ValueError(
                "dynamic_atoms is missing. Regenerate utm_scenario.json using "
                "the operation-composition/multispeed generator."
            )
        if not isinstance(atoms, list):
            continue

        for atom in atoms:
            if not isinstance(atom, Mapping):
                continue
            key = (int(atom["x21"]), int(atom["y21"]), int(atom["h12"]))
            intervals = atom.get("intervals", [])
            for interval in intervals:
                if not isinstance(interval, (list, tuple)) or len(interval) != 2:
                    continue
                start, end = int(interval[0]), int(interval[1])
                if end > start:
                    entries_by_atom[key].append((start, end, flight_id))

    return entries_by_atom, x_zoom, y_zoom, h_zoom


def write_csv(path: Path, rows: Iterable[Mapping[str, object]], fieldnames: Sequence[str]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def analyze(input_path: Path, output_dir: Path, dump_records: bool = True) -> Dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {input_path} ...")
    data = json.loads(input_path.read_text(encoding="utf-8"))

    flights = data.get("flights", [])
    if not isinstance(flights, list):
        flights = []
    flight_ids = [str(f["flight_id"]) for f in flights if isinstance(f, Mapping) and "flight_id" in f]
    flight_count = len(flight_ids)
    attribute_mask_bits = max(1, flight_count)

    entries_by_atom, zx0, zy0, zh0 = get_dynamic_entries(data)

    fixed_xyh_bits = bst_bits_3d(zx0, zy0, zh0)
    sampled_key_bits_per_record = fixed_xyh_bits + 32
    fixed_interval_key_bits_per_record = fixed_xyh_bits + 64

    print(f"Fixed working grid: X={zx0}, Y={zy0}, H={zh0}")
    print(f"Unique fixed XYH atoms: {len(entries_by_atom)}")
    print(f"Attribute mask model: {attribute_mask_bits} bits/record ({flight_count} flights)")

    # Full-record dumps are streamed to gzip so the analysis can retain the full
    # reproducible dataset without creating enormous plain CSV files.
    fixed_gz = None
    fixed_writer = None
    if dump_records:
        fixed_gz = gzip.open(output_dir / "fixed_interval_records.csv.gz", "wt", newline="", encoding="utf-8")
        fixed_writer = csv.DictWriter(
            fixed_gz,
            fieldnames=[
                "x_zoom", "x_index", "y_zoom", "y_index", "h_zoom", "h_index",
                "start_epoch_s", "end_epoch_s", "duration_s",
                "attribute_cardinality", "attribute_flights",
            ],
        )
        fixed_writer.writeheader()

    cells_by_attribute: Dict[Attribute, Set[Prefix]] = defaultdict(set)

    fixed_interval_count = 0
    fixed_interval_key_bits = 0
    fixed_interval_attribute_bits = 0
    sampled_cell_seconds = 0
    sampled_key_bits = 0
    sampled_attribute_bits = 0

    fixed_interval_by_cardinality: Counter[int] = Counter()
    fixed_dyadic_by_cardinality: Counter[int] = Counter()

    fixed_dyadic_count = 0
    fixed_dyadic_key_bits = 0
    fixed_dyadic_attribute_bits = 0

    print("Sweeping exact active-flight attributes and building dyadic T prefixes ...")
    for atom_idx, ((ix, iy, ih), entries) in enumerate(entries_by_atom.items(), 1):
        for start_s, end_s, attr in sweep_attribute_segments(entries):
            duration = end_s - start_s
            card = len(attr)

            fixed_interval_count += 1
            fixed_interval_key_bits += fixed_interval_key_bits_per_record
            fixed_interval_attribute_bits += attribute_mask_bits
            fixed_interval_by_cardinality[card] += 1

            sampled_cell_seconds += duration
            sampled_key_bits += duration * sampled_key_bits_per_record
            sampled_attribute_bits += duration * attribute_mask_bits

            if fixed_writer is not None:
                fixed_writer.writerow({
                    "x_zoom": zx0, "x_index": ix,
                    "y_zoom": zy0, "y_index": iy,
                    "h_zoom": zh0, "h_index": ih,
                    "start_epoch_s": start_s,
                    "end_epoch_s": end_s,
                    "duration_s": duration,
                    "attribute_cardinality": card,
                    "attribute_flights": "|".join(attr),
                })

            for zt, it in dyadic_time_cover(start_s, end_s):
                prefix: Prefix = (
                    (zx0, ix),
                    (zy0, iy),
                    (zh0, ih),
                    (zt, it),
                )
                if prefix not in cells_by_attribute[attr]:
                    cells_by_attribute[attr].add(prefix)
                    fixed_dyadic_count += 1
                    fixed_dyadic_key_bits += bst_bits_4d(zx0, zy0, zh0, zt)
                    fixed_dyadic_attribute_bits += attribute_mask_bits
                    fixed_dyadic_by_cardinality[card] += 1

        if atom_idx % 5000 == 0:
            print(f"  processed {atom_idx}/{len(entries_by_atom)} fixed atoms")

    if fixed_gz is not None:
        fixed_gz.close()

    print(
        f"Attribute groups: {len(cells_by_attribute)}; "
        f"fixed dyadic XYHT records: {fixed_dyadic_count}"
    )

    normalized_count = 0
    normalized_key_bits = 0
    normalized_attribute_bits = 0
    normalized_by_cardinality: Counter[int] = Counter()
    zoom_distribution: Counter[Tuple[int, int, int, int, int]] = Counter()
    merge_counts = Counter({"x": 0, "y": 0, "h": 0, "t": 0})
    coarsened_axis_records = Counter({"x": 0, "y": 0, "h": 0, "t": 0})
    anisotropic_xy_count = 0
    spatially_merged_count = 0

    normalized_gz = None
    normalized_writer = None
    if dump_records:
        normalized_gz = gzip.open(output_dir / "normalized_xyht_records.csv.gz", "wt", newline="", encoding="utf-8")
        normalized_writer = csv.DictWriter(
            normalized_gz,
            fieldnames=[
                "attribute_cardinality", "attribute_flights",
                "zx", "ix", "zy", "iy", "zh", "ih", "zt", "it",
                "time_block_s", "bst_id_bits",
            ],
        )
        normalized_writer.writeheader()

    print("Normalizing each exact attribute class in XYHT ...")
    attr_items = sorted(cells_by_attribute.items(), key=lambda kv: (len(kv[0]), kv[0]))
    for attr_idx, (attr, cells) in enumerate(attr_items, 1):
        normalized, merges = normalize_prefix_partition(cells)
        merge_counts.update(merges)
        card = len(attr)

        for prefix in normalized:
            (zx, ix), (zy, iy), (zh, ih), (zt, it) = prefix
            bits = bst_bits_4d(zx, zy, zh, zt)

            normalized_count += 1
            normalized_key_bits += bits
            normalized_attribute_bits += attribute_mask_bits
            normalized_by_cardinality[card] += 1
            zoom_distribution[(card, zx, zy, zh, zt)] += 1

            if zx < zx0:
                coarsened_axis_records["x"] += 1
            if zy < zy0:
                coarsened_axis_records["y"] += 1
            if zh < zh0:
                coarsened_axis_records["h"] += 1
            if zt < 32:
                coarsened_axis_records["t"] += 1
            if zx != zy:
                anisotropic_xy_count += 1
            if zx < zx0 or zy < zy0 or zh < zh0:
                spatially_merged_count += 1

            if normalized_writer is not None:
                normalized_writer.writerow({
                    "attribute_cardinality": card,
                    "attribute_flights": "|".join(attr),
                    "zx": zx, "ix": ix,
                    "zy": zy, "iy": iy,
                    "zh": zh, "ih": ih,
                    "zt": zt, "it": it,
                    "time_block_s": 1 << (T_BITS - zt),
                    "bst_id_bits": bits,
                })

        if attr_idx % 100 == 0:
            print(f"  normalized {attr_idx}/{len(attr_items)} attribute groups")

    if normalized_gz is not None:
        normalized_gz.close()

    # Consistency: each binary merge reduces the record count by exactly one.
    total_merges = sum(merge_counts.values())
    expected_normalized = fixed_dyadic_count - total_merges
    if expected_normalized != normalized_count:
        raise RuntimeError(
            f"normalization accounting mismatch: "
            f"{fixed_dyadic_count} - {total_merges} != {normalized_count}"
        )

    representations = [
        {
            "representation": "sampled_fixed_grid_1s",
            "records": sampled_cell_seconds,
            "key_bits": sampled_key_bits,
            "attribute_bits": sampled_attribute_bits,
            "total_bits": sampled_key_bits + sampled_attribute_bits,
            "definition": (
                f"one z{zx0}/z{zy0}/H{zh0} XYH atom + 32-bit timestamp "
                "for every active cell-second"
            ),
        },
        {
            "representation": "fixed_grid_exact_intervals",
            "records": fixed_interval_count,
            "key_bits": fixed_interval_key_bits,
            "attribute_bits": fixed_interval_attribute_bits,
            "total_bits": fixed_interval_key_bits + fixed_interval_attribute_bits,
            "definition": (
                f"fixed z{zx0}/z{zy0}/H{zh0} XYH atoms with maximal exact "
                "start/end uint32 intervals; no spatial merging"
            ),
        },
        {
            "representation": "fixed_grid_dyadic_T",
            "records": fixed_dyadic_count,
            "key_bits": fixed_dyadic_key_bits,
            "attribute_bits": fixed_dyadic_attribute_bits,
            "total_bits": fixed_dyadic_key_bits + fixed_dyadic_attribute_bits,
            "definition": (
                f"fixed z{zx0}/z{zy0}/H{zh0} XYH atoms with exact dyadic BST-ID "
                "T-prefix cover; no spatial merging"
            ),
        },
        {
            "representation": "BST_ID_normalized_XYHT",
            "records": normalized_count,
            "key_bits": normalized_key_bits,
            "attribute_bits": normalized_attribute_bits,
            "total_bits": normalized_key_bits + normalized_attribute_bits,
            "definition": (
                "exact-attribute XYHT prefix sets normalized by independent "
                "x/y/h/t binary sibling merging"
            ),
        },
    ]

    baseline_dyadic = next(r for r in representations if r["representation"] == "fixed_grid_dyadic_T")
    baseline_interval = next(r for r in representations if r["representation"] == "fixed_grid_exact_intervals")
    proposed = next(r for r in representations if r["representation"] == "BST_ID_normalized_XYHT")

    for row in representations:
        total = int(row["total_bits"])
        key = int(row["key_bits"])
        row["key_KiB"] = round(key / 8 / 1024, 6)
        row["total_KiB"] = round(total / 8 / 1024, 6)
        row["record_ratio_vs_fixed_dyadic"] = (
            round(int(row["records"]) / int(baseline_dyadic["records"]), 8)
            if int(baseline_dyadic["records"]) else None
        )
        row["total_bit_ratio_vs_fixed_dyadic"] = (
            round(total / int(baseline_dyadic["total_bits"]), 8)
            if int(baseline_dyadic["total_bits"]) else None
        )

    summary = {
        "input_file": str(input_path),
        "input_file_bytes": input_path.stat().st_size,
        "flight_count": flight_count,
        "fixed_grid": {"x_zoom": zx0, "y_zoom": zy0, "h_zoom": zh0},
        "attribute_model": {
            "semantic_rule": "cells/prefixes merge only when the exact attribute value is identical",
            "scenario_attribute": "exact active reservation flight set",
            "comparison_encoding": f"{attribute_mask_bits}-bit flight-set bitmask per record",
            "note": (
                "The flight-set attribute is used only as a controlled scenario encoding. "
                "Application-specific attribute dictionaries/metadata are excluded."
            ),
        },
        "representations": representations,
        "proposed_vs_fixed_dyadic": {
            "record_reduction_fraction": (
                1.0 - normalized_count / fixed_dyadic_count if fixed_dyadic_count else 0.0
            ),
            "record_compression_factor": (
                fixed_dyadic_count / normalized_count if normalized_count else None
            ),
            "key_bit_reduction_fraction": (
                1.0 - normalized_key_bits / fixed_dyadic_key_bits if fixed_dyadic_key_bits else 0.0
            ),
            "key_bit_compression_factor": (
                fixed_dyadic_key_bits / normalized_key_bits if normalized_key_bits else None
            ),
            "total_bit_reduction_fraction_with_attribute_mask": (
                1.0 - (normalized_key_bits + normalized_attribute_bits)
                / (fixed_dyadic_key_bits + fixed_dyadic_attribute_bits)
                if (fixed_dyadic_key_bits + fixed_dyadic_attribute_bits) else 0.0
            ),
            "total_bit_compression_factor_with_attribute_mask": (
                (fixed_dyadic_key_bits + fixed_dyadic_attribute_bits)
                / (normalized_key_bits + normalized_attribute_bits)
                if (normalized_key_bits + normalized_attribute_bits) else None
            ),
        },
        "proposed_vs_fixed_exact_intervals": {
            "record_ratio": normalized_count / fixed_interval_count if fixed_interval_count else None,
            "key_bit_ratio": normalized_key_bits / fixed_interval_key_bits if fixed_interval_key_bits else None,
            "note": (
                "The exact-interval baseline is intentionally strong: it permits arbitrary "
                "start/end intervals, whereas BST-ID uses dyadic T prefixes."
            ),
        },
        "normalization": {
            "policy": (
                "deterministic x->y->h->t binary sibling merging; restart from x "
                "after each merge layer; no claim of global minimum cover"
            ),
            "binary_merge_count_total": total_merges,
            "binary_merge_count_by_axis": dict(merge_counts),
            "final_spatially_merged_records": spatially_merged_count,
            "final_spatially_merged_fraction": (
                spatially_merged_count / normalized_count if normalized_count else 0.0
            ),
            "final_anisotropic_xy_records": anisotropic_xy_count,
            "final_anisotropic_xy_fraction": (
                anisotropic_xy_count / normalized_count if normalized_count else 0.0
            ),
            "final_records_coarsened_by_axis": dict(coarsened_axis_records),
        },
    }

    # Main summary JSON.
    (output_dir / "compression_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    # Human-readable report for quick inspection and paper drafting.
    proposed_metrics = summary["proposed_vs_fixed_dyadic"]
    norm_metrics = summary["normalization"]
    report = f"""# BST-ID whole-scenario compression analysis

## Evaluation scope

- Flights: {flight_count}
- Fixed working grid: X z{zx0}, Y z{zy0}, H z{zh0}
- Semantic attribute: exact active reservation-flight set
- Merge rule: prefixes are mergeable only when their semantic attribute is identical.
- Attribute comparison model: {attribute_mask_bits}-bit flight-set mask per record.

The main comparison isolates the effect of BST-ID prefix normalization:

1. **Fixed grid + dyadic T**: XYH remains fixed at z{zx0}/z{zy0}/H{zh0}; exact time intervals are represented by BST-ID-compatible dyadic T prefixes.
2. **BST-ID normalized XYHT**: the same semantic data are normalized by independent binary sibling merging along X, Y, H, and T.

The arbitrary exact-interval baseline is also reported as a strong practical reference. It can encode any start/end pair directly, whereas BST-ID T prefixes are dyadic.

## Main result

- Fixed-grid dyadic records: **{fixed_dyadic_count:,}**
- Normalized BST-ID records: **{normalized_count:,}**
- Record reduction: **{proposed_metrics['record_reduction_fraction'] * 100:.2f}%**
- Record compression factor: **{proposed_metrics['record_compression_factor']:.3f}x**
- Key-bit reduction: **{proposed_metrics['key_bit_reduction_fraction'] * 100:.2f}%**
- Key-bit compression factor: **{proposed_metrics['key_bit_compression_factor']:.3f}x**
- Total-bit reduction including the controlled attribute mask: **{proposed_metrics['total_bit_reduction_fraction_with_attribute_mask'] * 100:.2f}%**

## Why records were reduced

Binary sibling merges:

- X-axis: **{merge_counts['x']:,}**
- Y-axis: **{merge_counts['y']:,}**
- H-axis: **{merge_counts['h']:,}**
- T-axis: **{merge_counts['t']:,}**

Final normalized records with at least one spatial axis coarsened:
**{spatially_merged_count:,}/{normalized_count:,} ({norm_metrics['final_spatially_merged_fraction'] * 100:.2f}%)**

Final records with different X/Y zooms (direct evidence of rectangular, non-square horizontal prefixes):
**{anisotropic_xy_count:,}/{normalized_count:,} ({norm_metrics['final_anisotropic_xy_fraction'] * 100:.2f}%)**

## Output files

- `compression_summary.json`: complete machine-readable summary
- `compression_summary.csv`: representation-level comparison
- `axis_merge_summary.csv`: merge counts by X/Y/H/T
- `zoom_distribution.csv`: final zoom-vector distribution
- `attribute_cardinality_distribution.csv`: ordinary vs overlap-region statistics
- `flight_summary.csv`: per-flight context including 5/8/10 m/s cruise speed
- `fixed_interval_records.csv.gz`: full fixed-grid exact-interval records
- `normalized_xyht_records.csv.gz`: full normalized XYHT prefix records

## Interpretation note

This analysis does **not** claim a globally minimum rectangle cover. It applies a deterministic X->Y->H->T binary sibling-normalization policy while preserving the exact semantic attribute. The comparison therefore evaluates the actual hierarchical representation produced by the stated policy, not an unconstrained optimal packing.
"""
    (output_dir / "compression_report.md").write_text(report, encoding="utf-8")

    # Compact representation comparison CSV.
    write_csv(
        output_dir / "compression_summary.csv",
        representations,
        [
            "representation", "records", "key_bits", "attribute_bits", "total_bits",
            "key_KiB", "total_KiB",
            "record_ratio_vs_fixed_dyadic", "total_bit_ratio_vs_fixed_dyadic",
            "definition",
        ],
    )

    # Axis merge counts.
    write_csv(
        output_dir / "axis_merge_summary.csv",
        [
            {
                "axis": axis,
                "binary_sibling_merges": merge_counts[axis],
                "fraction_of_all_merges": (
                    merge_counts[axis] / total_merges if total_merges else 0.0
                ),
            }
            for axis in ("x", "y", "h", "t")
        ],
        ["axis", "binary_sibling_merges", "fraction_of_all_merges"],
    )

    # Zoom-vector distribution, retaining attribute cardinality so overlap
    # regions can be studied without tying the analysis to a particular color.
    zoom_rows = []
    for (card, zx, zy, zh, zt), count in sorted(zoom_distribution.items()):
        zoom_rows.append({
            "attribute_cardinality": card,
            "zx": zx, "zy": zy, "zh": zh, "zt": zt,
            "count": count,
            "bst_bits_per_record": bst_bits_4d(zx, zy, zh, zt),
            "time_block_s": 1 << (T_BITS - zt),
            "x_coarsened_bits": zx0 - zx,
            "y_coarsened_bits": zy0 - zy,
            "h_coarsened_bits": zh0 - zh,
            "anisotropic_xy": int(zx != zy),
        })
    write_csv(
        output_dir / "zoom_distribution.csv",
        zoom_rows,
        [
            "attribute_cardinality", "zx", "zy", "zh", "zt", "count",
            "bst_bits_per_record", "time_block_s",
            "x_coarsened_bits", "y_coarsened_bits", "h_coarsened_bits",
            "anisotropic_xy",
        ],
    )

    # Attribute cardinality distribution (1 = ordinary reservation, 2+ = overlap).
    cards = sorted(
        set(fixed_interval_by_cardinality)
        | set(fixed_dyadic_by_cardinality)
        | set(normalized_by_cardinality)
    )
    write_csv(
        output_dir / "attribute_cardinality_distribution.csv",
        [
            {
                "attribute_cardinality": c,
                "fixed_interval_records": fixed_interval_by_cardinality[c],
                "fixed_dyadic_records": fixed_dyadic_by_cardinality[c],
                "normalized_records": normalized_by_cardinality[c],
                "normalized_vs_fixed_dyadic_record_ratio": (
                    normalized_by_cardinality[c] / fixed_dyadic_by_cardinality[c]
                    if fixed_dyadic_by_cardinality[c] else None
                ),
            }
            for c in cards
        ],
        [
            "attribute_cardinality",
            "fixed_interval_records",
            "fixed_dyadic_records",
            "normalized_records",
            "normalized_vs_fixed_dyadic_record_ratio",
        ],
    )

    # Per-flight context, including the new 5/8/10 m/s assignment.
    flight_rows = []
    for flight in flights:
        if not isinstance(flight, Mapping):
            continue
        conf = flight.get("conformance", {})
        if not isinstance(conf, Mapping):
            conf = {}
        route = flight.get("route", {})
        if not isinstance(route, Mapping):
            route = {}
        reservation = flight.get("reservation", {})
        if not isinstance(reservation, Mapping):
            reservation = {}
        applicable = int(conf.get("applicable", conf.get("applicable_samples", 0)) or 0)
        inside = int(conf.get("inside", conf.get("inside_samples", 0)) or 0)
        flight_rows.append({
            "flight_id": flight.get("flight_id"),
            "origin": flight.get("origin_name"),
            "destination": flight.get("destination_name"),
            "target_agl_m": flight.get("target_agl_m"),
            "planned_cruise_speed_mps": flight.get(
                "planned_cruise_speed_mps",
                route.get("cruise_speed_mps"),
            ),
            "distance_m": route.get("distance_m"),
            "nominal_duration_s": route.get("nominal_duration_s"),
            "cruise_duration_s": route.get("cruise_duration_s"),
            "dynamic_atom_count": reservation.get("dynamic_atom_count"),
            "display_cell_count": reservation.get("display_cell_count"),
            "conformance_applicable_samples": applicable,
            "conformance_inside_samples": inside,
            "conformance_inside_ratio": (inside / applicable if applicable else None),
        })
    write_csv(
        output_dir / "flight_summary.csv",
        flight_rows,
        [
            "flight_id", "origin", "destination", "target_agl_m",
            "planned_cruise_speed_mps", "distance_m",
            "nominal_duration_s", "cruise_duration_s",
            "dynamic_atom_count", "display_cell_count",
            "conformance_applicable_samples", "conformance_inside_samples",
            "conformance_inside_ratio",
        ],
    )

    print()
    print("Compression analysis complete")
    print("-----------------------------")
    print(f"Fixed exact interval records : {fixed_interval_count:,}")
    print(f"Fixed dyadic XYHT records    : {fixed_dyadic_count:,}")
    print(f"Normalized BST-ID records    : {normalized_count:,}")
    if normalized_count:
        print(f"Record compression factor    : {fixed_dyadic_count / normalized_count:.3f}x")
    if normalized_key_bits:
        print(f"Key-bit compression factor   : {fixed_dyadic_key_bits / normalized_key_bits:.3f}x")
    if normalized_count:
        print(
            f"Anisotropic XY prefixes      : {anisotropic_xy_count:,}/"
            f"{normalized_count:,} "
            f"({anisotropic_xy_count / normalized_count * 100:.2f}%)"
        )
    print(
        "Sibling merges x/y/h/t     : "
        f"{merge_counts['x']:,} / {merge_counts['y']:,} / "
        f"{merge_counts['h']:,} / {merge_counts['t']:,}"
    )
    print(f"Output directory             : {output_dir}")
    if dump_records:
        print("Full records                 : *.csv.gz included")
    print()

    return summary


def main() -> None:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=(
            "Analyze whole-scenario BST-ID reservation compression with exact "
            "attribute preservation."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=here / "utm_scenario.json",
        help="Path to utm_scenario.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=here / "compression_analysis",
        help="Directory for JSON/CSV outputs",
    )
    parser.add_argument(
        "--no-record-dumps",
        action="store_true",
        help="Skip the two full *.csv.gz record dumps",
    )
    args = parser.parse_args()

    analyze(
        input_path=args.input.resolve(),
        output_dir=args.output_dir.resolve(),
        dump_records=not args.no_record_dumps,
    )


if __name__ == "__main__":
    main()
