#!/usr/bin/env python3
"""
check_tsuruga_1m_pipeline.py

Diagnostic checker for the Tsuruga Castle BST-ID-native pipeline.

It reconstructs the full pipeline from the aligned point cloud:

  aligned point cloud
    -> raw fine BST-ID shell
    -> 3x3x3 closing
    -> directional fill
    -> solid fine BST-ID object
    -> normalized-prefix expansion

and checks where geometry is being changed.

Outputs:
  <prefix>_report.json
  <prefix>_classified.ply
  <prefix>_raw_shell.ply
  <prefix>_closed_added.ply
  <prefix>_fill_added.ply
  <prefix>_fill_exposed_boundary.ply
  <prefix>_normalize_missing.ply
  <prefix>_normalize_extra.ply
  <prefix>_fill_exposed_boundary.csv

The most useful quantity is:

  final solid boundary cells NOT present in the closed shell

Those cells were created by the fill rule and became part of the exterior
surface. If that number is large, the fill rule is shaping the outside of the
object instead of only filling an enclosed interior.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import open3d as o3d
from pyproj import CRS, Transformer
from scipy import ndimage


WEB_MERCATOR_LAT_MAX = 85.0511287798066
H_MIN = -16383.0
H_MAX = 16384.0
H_SPAN = H_MAX - H_MIN


def load_points(path: Path) -> np.ndarray:
    pcd = o3d.io.read_point_cloud(str(path))
    if pcd.is_empty():
        raise RuntimeError(f"Could not read point cloud: {path}")
    return np.asarray(pcd.points).copy()


def make_local_to_wgs84(anchor_lat: float, anchor_lon: float):
    local = CRS.from_proj4(
        f"+proj=aeqd +lat_0={anchor_lat:.12f} "
        f"+lon_0={anchor_lon:.12f} +datum=WGS84 +units=m +no_defs"
    )
    return Transformer.from_crs(local, CRS.from_epsg(4326), always_xy=True)


def make_wgs84_to_local(anchor_lat: float, anchor_lon: float):
    local = CRS.from_proj4(
        f"+proj=aeqd +lat_0={anchor_lat:.12f} "
        f"+lon_0={anchor_lon:.12f} +datum=WGS84 +units=m +no_defs"
    )
    return Transformer.from_crs(CRS.from_epsg(4326), local, always_xy=True)


def lonlat_to_prefix(lon_deg, lat_deg, zx: int, zy: int):
    lon = np.asarray(lon_deg, dtype=np.float64)
    lat = np.asarray(lat_deg, dtype=np.float64)

    ix = np.floor(((lon + 180.0) / 360.0) * (1 << zx)).astype(np.int64)
    ix = np.clip(ix, 0, (1 << zx) - 1)

    lat = np.clip(lat, -WEB_MERCATOR_LAT_MAX, WEB_MERCATOR_LAT_MAX)
    lat_rad = np.radians(lat)
    yn = (1.0 - np.arcsinh(np.tan(lat_rad)) / math.pi) / 2.0

    iy = np.floor(yn * (1 << zy)).astype(np.int64)
    iy = np.clip(iy, 0, (1 << zy) - 1)
    return ix, iy


def altitude_to_prefix(h_m, zh: int):
    h = np.asarray(h_m, dtype=np.float64)
    u = (h - H_MIN) / H_SPAN
    u = np.clip(u, 0.0, np.nextafter(1.0, 0.0))
    return np.floor(u * (1 << zh)).astype(np.int64)


def direct_quantize(
    points_enu: np.ndarray,
    anchor_lat: float,
    anchor_lon: float,
    height_origin: float,
    zx: int,
    zy: int,
    zh: int,
):
    to_wgs84 = make_local_to_wgs84(anchor_lat, anchor_lon)
    lon, lat = to_wgs84.transform(points_enu[:, 0], points_enu[:, 1])

    ix, iy = lonlat_to_prefix(lon, lat, zx, zy)
    ih = altitude_to_prefix(height_origin + points_enu[:, 2], zh)

    return np.unique(np.column_stack((ix, iy, ih)), axis=0)


def indices_to_dense(indices: np.ndarray, margin: int):
    mins = indices.min(axis=0).astype(np.int64) - margin
    maxs = indices.max(axis=0).astype(np.int64) + margin
    shape = (maxs - mins + 1).astype(int)

    grid = np.zeros(tuple(shape), dtype=bool)
    local = indices - mins
    grid[local[:, 0], local[:, 1], local[:, 2]] = True
    return grid, mins


def dense_to_indices(mask: np.ndarray, mins: np.ndarray):
    return np.argwhere(mask).astype(np.int64) + mins


def directional_fill(shell: np.ndarray, require_ns: bool):
    west = np.zeros_like(shell)
    c = np.maximum.accumulate(shell, axis=0)
    west[1:, :, :] = c[:-1, :, :]

    east = np.zeros_like(shell)
    c = np.maximum.accumulate(shell[::-1, :, :], axis=0)[::-1, :, :]
    east[:-1, :, :] = c[1:, :, :]

    up = np.zeros_like(shell)
    c = np.maximum.accumulate(shell[:, :, ::-1], axis=2)[:, :, ::-1]
    up[:, :, :-1] = c[:, :, 1:]

    interior = (~shell) & west & east & up

    if require_ns:
        side_a = np.zeros_like(shell)
        c = np.maximum.accumulate(shell, axis=1)
        side_a[:, 1:, :] = c[:, :-1, :]

        side_b = np.zeros_like(shell)
        c = np.maximum.accumulate(shell[:, ::-1, :], axis=1)[:, ::-1, :]
        side_b[:, :-1, :] = c[:, 1:, :]

        interior &= side_a & side_b

    return shell | interior, interior


def load_fine_from_cesium(path: Path):
    d = json.loads(path.read_text(encoding="utf-8"))
    return {
        (int(r["ix"]), int(r["iy"]), int(r["ih"]))
        for r in d["bitmap"]
    }


def load_and_expand_normalized(path: Path, zx: int, zy: int, zh: int):
    expanded = set()
    rows = []

    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            row = {k: int(v) for k, v in r.items()}
            rows.append(row)

    for r in rows:
        dx = zx - r["zx"]
        dy = zy - r["zy"]
        dh = zh - r["zh"]

        x0 = r["ix"] << dx
        y0 = r["iy"] << dy
        h0 = r["ih"] << dh

        for ix in range(x0, x0 + (1 << dx)):
            for iy in range(y0, y0 + (1 << dy)):
                for ih in range(h0, h0 + (1 << dh)):
                    expanded.add((ix, iy, ih))

    return rows, expanded


def y_to_lat(y: float, z: int):
    n = float(1 << z)
    v = math.pi * (1.0 - 2.0 * y / n)
    return math.degrees(math.atan(math.sinh(v)))


def fine_indices_to_local_centers(
    indices,
    zx: int,
    zy: int,
    zh: int,
    anchor_lat: float,
    anchor_lon: float,
    height_origin: float,
):
    if not indices:
        return np.empty((0, 3), dtype=np.float64)

    a = np.asarray(list(indices), dtype=np.int64)

    lon = ((a[:, 0].astype(np.float64) + 0.5) / (1 << zx)) * 360.0 - 180.0
    lat = np.array(
        [y_to_lat(float(v) + 0.5, zy) for v in a[:, 1]],
        dtype=np.float64,
    )

    to_local = make_wgs84_to_local(anchor_lat, anchor_lon)
    east, north = to_local.transform(lon, lat)

    h_width = H_SPAN / float(1 << zh)
    h_abs = H_MIN + (a[:, 2].astype(np.float64) + 0.5) * h_width
    up = h_abs - height_origin

    return np.column_stack((east, north, up))


def save_point_set(
    path: Path,
    cells,
    *,
    zx: int,
    zy: int,
    zh: int,
    anchor_lat: float,
    anchor_lon: float,
    height_origin: float,
    color=None,
):
    centers = fine_indices_to_local_centers(
        cells, zx, zy, zh, anchor_lat, anchor_lon, height_origin
    )

    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(centers)

    if color is not None and len(centers):
        colors = np.tile(np.asarray(color, dtype=np.float64), (len(centers), 1))
        pcd.colors = o3d.utility.Vector3dVector(colors)

    o3d.io.write_point_cloud(str(path), pcd, write_ascii=False)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--aligned-ply",
        type=Path,
        default=Path("tsuruga_aligned.ply"),
    )
    ap.add_argument(
        "--stats",
        type=Path,
        default=Path("tsuruga_native_1m_stats.json"),
    )
    ap.add_argument(
        "--cesium",
        type=Path,
        default=Path("tsuruga_native_1m_cesium.json"),
    )
    ap.add_argument(
        "--normalized",
        type=Path,
        default=Path("tsuruga_native_1m_normalized.csv"),
    )
    ap.add_argument(
        "--out-prefix",
        default="tsuruga_diag",
    )

    args = ap.parse_args()

    stats = json.loads(args.stats.read_text(encoding="utf-8"))

    zx = int(stats["fine_grid"]["zx"])
    zy = int(stats["fine_grid"]["zy"])
    zh = int(stats["fine_grid"]["zh"])

    anchor_lat = float(stats["georeference"]["anchor_latitude_deg"])
    anchor_lon = float(stats["georeference"]["anchor_longitude_deg"])
    height_origin = float(stats["georeference"]["height_origin_m"])

    min_u = float(stats["input"]["min_u_m"])
    margin = 2
    close_iterations = int(stats["fine_grid"]["close_iterations"])
    require_ns = bool(stats["fine_grid"].get("require_ns", False))

    points = load_points(args.aligned_ply)
    input_point_count = len(points)

    points = points[points[:, 2] >= min_u]

    raw_idx = direct_quantize(
        points,
        anchor_lat,
        anchor_lon,
        height_origin,
        zx,
        zy,
        zh,
    )

    raw_grid, dense_origin = indices_to_dense(raw_idx, margin)

    closed_grid = ndimage.binary_closing(
        raw_grid,
        structure=np.ones((3, 3, 3), dtype=bool),
        iterations=close_iterations,
    )

    solid_grid, interior_grid = directional_fill(
        closed_grid,
        require_ns=require_ns,
    )

    raw = set(map(tuple, raw_idx))
    closed = set(map(tuple, dense_to_indices(closed_grid, dense_origin)))
    solid = set(map(tuple, dense_to_indices(solid_grid, dense_origin)))
    interior = set(map(tuple, dense_to_indices(interior_grid, dense_origin)))

    current_fine = load_fine_from_cesium(args.cesium)
    normalized_rows, normalized_expanded = load_and_expand_normalized(
        args.normalized, zx, zy, zh
    )

    # Final 6-neighbour exterior boundary of the solid.
    structure6 = ndimage.generate_binary_structure(rank=3, connectivity=1)
    eroded = ndimage.binary_erosion(
        solid_grid,
        structure=structure6,
        border_value=0,
    )
    boundary_grid = solid_grid & (~eroded)
    boundary = set(
        map(tuple, dense_to_indices(boundary_grid, dense_origin))
    )

    # These cells are especially interesting:
    # they are on the FINAL exterior surface but were NOT in the closed shell.
    # Therefore they were created by the fill rule.
    fill_exposed_boundary = boundary - closed

    closed_added = closed - raw
    fill_added = solid - closed

    normalize_missing = current_fine - normalized_expanded
    normalize_extra = normalized_expanded - current_fine

    pipeline_missing = solid - current_fine
    pipeline_extra = current_fine - solid

    # Exposed-face classification for fill-created boundary cells.
    dirs = {
        "+x": (1, 0, 0),
        "-x": (-1, 0, 0),
        "+y": (0, 1, 0),
        "-y": (0, -1, 0),
        "+h": (0, 0, 1),
        "-h": (0, 0, -1),
    }

    exposed_face_counts = Counter()
    exposed_rows = []

    for p in sorted(fill_exposed_boundary):
        faces = []

        for name, d in dirs.items():
            q = (p[0] + d[0], p[1] + d[1], p[2] + d[2])
            if q not in solid:
                exposed_face_counts[name] += 1
                faces.append(name)

        exposed_rows.append(
            {
                "ix": p[0],
                "iy": p[1],
                "ih": p[2],
                "exposed_faces": " ".join(faces),
            }
        )

    report = {
        "input": {
            "aligned_ply": str(args.aligned_ply),
            "input_point_count": input_point_count,
            "expected_input_point_count_from_stats":
                int(stats["input"]["input_points"]),
            "points_after_u_crop": len(points),
            "expected_points_after_crop_from_stats":
                int(stats["input"]["points_after_crop"]),
        },
        "pipeline": {
            "raw_shell_cells": len(raw),
            "expected_raw_shell_cells":
                int(stats["fine_grid"]["raw_shell_cells"]),
            "closed_shell_cells": len(closed),
            "expected_closed_shell_cells":
                int(stats["fine_grid"]["closed_shell_cells"]),
            "closed_added_cells": len(closed_added),
            "closed_removed_raw_cells": len(raw - closed),
            "fill_added_cells": len(fill_added),
            "solid_cells": len(solid),
            "expected_solid_cells":
                int(stats["fine_grid"]["solid_fine_cells"]),
            "current_cesium_fine_cells": len(current_fine),
            "reconstructed_vs_current_missing": len(pipeline_missing),
            "reconstructed_vs_current_extra": len(pipeline_extra),
        },
        "normalize": {
            "normalized_prefix_count": len(normalized_rows),
            "expanded_fine_cells": len(normalized_expanded),
            "missing_fine_cells": len(normalize_missing),
            "extra_fine_cells": len(normalize_extra),
            "exact": (
                len(normalize_missing) == 0
                and len(normalize_extra) == 0
            ),
        },
        "boundary_diagnostic": {
            "final_solid_boundary_cells": len(boundary),
            "boundary_cells_already_in_raw_shell": len(boundary & raw),
            "boundary_cells_in_closed_shell": len(boundary & closed),
            "boundary_cells_created_only_by_fill":
                len(fill_exposed_boundary),
            "fill_created_boundary_fraction":
                (
                    len(fill_exposed_boundary) / len(boundary)
                    if boundary else 0.0
                ),
            "fill_created_exposed_face_counts":
                dict(exposed_face_counts),
        },
    }

    prefix = Path(args.out_prefix)

    Path(f"{prefix}_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    with Path(f"{prefix}_fill_exposed_boundary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["ix", "iy", "ih", "exposed_faces"],
        )
        writer.writeheader()
        writer.writerows(exposed_rows)

    # Individual PLY layers.
    save_point_set(
        Path(f"{prefix}_raw_shell.ply"),
        raw,
        zx=zx, zy=zy, zh=zh,
        anchor_lat=anchor_lat, anchor_lon=anchor_lon,
        height_origin=height_origin,
        color=[0.65, 0.65, 0.65],
    )

    save_point_set(
        Path(f"{prefix}_closed_added.ply"),
        closed_added,
        zx=zx, zy=zy, zh=zh,
        anchor_lat=anchor_lat, anchor_lon=anchor_lon,
        height_origin=height_origin,
        color=[1.0, 0.65, 0.0],
    )

    save_point_set(
        Path(f"{prefix}_fill_added.ply"),
        fill_added,
        zx=zx, zy=zy, zh=zh,
        anchor_lat=anchor_lat, anchor_lon=anchor_lon,
        height_origin=height_origin,
        color=[0.15, 0.45, 1.0],
    )

    save_point_set(
        Path(f"{prefix}_fill_exposed_boundary.ply"),
        fill_exposed_boundary,
        zx=zx, zy=zy, zh=zh,
        anchor_lat=anchor_lat, anchor_lon=anchor_lon,
        height_origin=height_origin,
        color=[1.0, 0.0, 0.0],
    )

    save_point_set(
        Path(f"{prefix}_normalize_missing.ply"),
        normalize_missing,
        zx=zx, zy=zy, zh=zh,
        anchor_lat=anchor_lat, anchor_lon=anchor_lon,
        height_origin=height_origin,
        color=[1.0, 0.0, 1.0],
    )

    save_point_set(
        Path(f"{prefix}_normalize_extra.ply"),
        normalize_extra,
        zx=zx, zy=zy, zh=zh,
        anchor_lat=anchor_lat, anchor_lon=anchor_lon,
        height_origin=height_origin,
        color=[0.0, 1.0, 1.0],
    )

    # One combined, color-coded diagnostic point cloud.
    #
    # Priority:
    #   red    = fill-created FINAL exterior boundary
    #   blue   = other fill-added cells
    #   orange = closing-added cells
    #   gray   = original raw shell
    categories = [
        (raw, [0.62, 0.62, 0.62]),
        (closed_added, [1.0, 0.62, 0.0]),
        (fill_added - fill_exposed_boundary, [0.15, 0.45, 1.0]),
        (fill_exposed_boundary, [1.0, 0.0, 0.0]),
    ]

    all_points = []
    all_colors = []

    for cells, color in categories:
        centers = fine_indices_to_local_centers(
            cells, zx, zy, zh, anchor_lat, anchor_lon, height_origin
        )
        if len(centers):
            all_points.append(centers)
            all_colors.append(
                np.tile(np.asarray(color), (len(centers), 1))
            )

    combined = o3d.geometry.PointCloud()

    if all_points:
        combined.points = o3d.utility.Vector3dVector(
            np.vstack(all_points)
        )
        combined.colors = o3d.utility.Vector3dVector(
            np.vstack(all_colors)
        )

    o3d.io.write_point_cloud(
        f"{prefix}_classified.ply",
        combined,
        write_ascii=False,
    )

    print()
    print("Tsuruga BST-ID pipeline diagnostic")
    print("=================================")
    print(
        f"Input points: {input_point_count} "
        f"(stats: {stats['input']['input_points']})"
    )
    print(
        f"Raw shell:    {len(raw)} "
        f"(stats: {stats['fine_grid']['raw_shell_cells']})"
    )
    print(
        f"Closed shell: {len(closed)} "
        f"(stats: {stats['fine_grid']['closed_shell_cells']})"
    )
    print(
        f"Solid:        {len(solid)} "
        f"(stats: {stats['fine_grid']['solid_fine_cells']})"
    )
    print()
    print(
        "Reconstructed solid vs saved fine occupancy: "
        f"missing={len(pipeline_missing)}, extra={len(pipeline_extra)}"
    )
    print(
        "Normalized expansion vs fine occupancy: "
        f"missing={len(normalize_missing)}, extra={len(normalize_extra)}"
    )
    print()
    print(f"Final solid boundary cells: {len(boundary)}")
    print(
        "Boundary created ONLY by directional fill: "
        f"{len(fill_exposed_boundary)} "
        f"({100.0 * report['boundary_diagnostic']['fill_created_boundary_fraction']:.2f}%)"
    )
    print(
        "Exposed faces of fill-created boundary:",
        dict(exposed_face_counts),
    )
    print()
    print("Most useful visualization:")
    print(f"  {prefix}_classified.ply")
    print()
    print("Color legend:")
    print("  gray   raw shell")
    print("  orange added by closing")
    print("  blue   added by directional fill, not exterior")
    print("  red    added by directional fill AND exposed on final exterior")
    print()
    print("If normalize_missing/extra are both 0, Normalize is not deleting geometry.")


if __name__ == "__main__":
    main()
