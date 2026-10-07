#!/usr/bin/env python3
"""
tsuruga_ply_to_bstid_native.py

BST-ID-native Tsuruga Castle object generation.

Key change from the earlier local-bitmap workflow:
    The aligned point cloud is FIRST quantized directly into the declared
    fine BST-ID XYH grid.  Surface closing, interior filling, and object
    construction are then performed in BST-ID integer-index space.

Pipeline
--------
aligned PLY (+X East, +Y North, +Z Up)
    -> crop local U >= min_u
    -> direct WGS84/BST-ID fine XYH quantization
    -> fine BST-ID shell
    -> morphology in (ix, iy, ih) index space
    -> directional or exterior-flood interior filling
    -> fine BST-ID object
    -> indexed/batch Normalize
    -> maximum inscribed BST-ID cell
    -> minimum enclosing BST-ID cell / valid cover
    -> Cesium JSON export

This avoids the grid-mismatch artifact caused by:
    local 1/2 m voxel centers -> independently aligned global BST-ID cells.

Run from the BST-ID repository root (or from journal/ if the package is
installed in the active venv).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import open3d as o3d
from pyproj import CRS, Transformer
from scipy import ndimage


try:
    from bst_id.region_algebra import Cell, normalize, region_atoms
except ImportError:
    # Typical layout when this script is placed in repository/journal/.
    repo_root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(repo_root))
    from bst_id.region_algebra import Cell, normalize, region_atoms

try:
    from bst_id.region_index import normalize_indexed_batch
    HAVE_INDEXED_BATCH = True
except ImportError:
    normalize_indexed_batch = None
    HAVE_INDEXED_BATCH = False


WEB_MERCATOR_LAT_MAX = 85.0511287798066
H_MIN = -16383.0
H_MAX = 16384.0
H_SPAN = H_MAX - H_MIN


# ---------------------------------------------------------------------------
# Coordinate / quantization utilities
# ---------------------------------------------------------------------------

def load_transform_metadata(path: Path) -> dict:
    d = json.loads(path.read_text(encoding="utf-8"))
    g = d["geographic_anchor"]
    return {
        "lat": float(g["latitude_deg"]),
        "lon": float(g["longitude_deg"]),
        "absolute_alt": g.get("absolute_altitude_m"),
        "raw": d,
    }


def make_local_to_wgs84(anchor_lat: float, anchor_lon: float):
    local_crs = CRS.from_proj4(
        f"+proj=aeqd +lat_0={anchor_lat:.12f} "
        f"+lon_0={anchor_lon:.12f} +datum=WGS84 +units=m +no_defs"
    )
    return Transformer.from_crs(
        local_crs,
        CRS.from_epsg(4326),
        always_xy=True,
    )


def make_wgs84_to_local(anchor_lat: float, anchor_lon: float):
    local_crs = CRS.from_proj4(
        f"+proj=aeqd +lat_0={anchor_lat:.12f} "
        f"+lon_0={anchor_lon:.12f} +datum=WGS84 +units=m +no_defs"
    )
    return Transformer.from_crs(
        CRS.from_epsg(4326),
        local_crs,
        always_xy=True,
    )


def lonlat_to_prefix(lon_deg, lat_deg, zx: int, zy: int):
    lon_deg = np.asarray(lon_deg, dtype=np.float64)
    lat_deg = np.asarray(lat_deg, dtype=np.float64)

    lon = ((lon_deg + 180.0) % 360.0) - 180.0
    lat = np.clip(lat_deg, -WEB_MERCATOR_LAT_MAX, WEB_MERCATOR_LAT_MAX)

    nx = float(1 << zx)
    ix = np.floor(((lon + 180.0) / 360.0) * nx).astype(np.int64)
    ix = np.clip(ix, 0, (1 << zx) - 1)

    ny = float(1 << zy)
    lat_rad = np.radians(lat)
    ynorm = (1.0 - np.arcsinh(np.tan(lat_rad)) / math.pi) / 2.0
    iy = np.floor(ynorm * ny).astype(np.int64)
    iy = np.clip(iy, 0, (1 << zy) - 1)

    return ix, iy


def altitude_to_prefix(h_m, zh: int):
    h_m = np.asarray(h_m, dtype=np.float64)
    n = float(1 << zh)
    u = (h_m - H_MIN) / H_SPAN
    u = np.clip(u, 0.0, np.nextafter(1.0, 0.0))
    return np.floor(u * n).astype(np.int64)


def height_boundary_near_zero(zh: int) -> float:
    """Choose a BST-ID H boundary to coincide with local U=0."""
    width = H_SPAN / float(1 << zh)
    k = round((0.0 - H_MIN) / width)
    k = max(0, min(1 << zh, k))
    return H_MIN + k * width


def resolve_height_origin(meta: dict, zh: int, user_anchor_alt):
    """
    Returns the height offset added to local U before H quantization.

    If a real absolute altitude is supplied, use it.
    Otherwise align local U=0 to the nearest BST-ID H cell boundary.  This
    prevents a fine H cell from straddling the chosen local ground plane.
    """
    if user_anchor_alt is not None:
        return float(user_anchor_alt), "command-line absolute altitude", False

    if meta["absolute_alt"] is not None:
        return float(meta["absolute_alt"]), "transform JSON absolute altitude", False

    offset = height_boundary_near_zero(zh)
    return (
        float(offset),
        "local relative height; U=0 aligned to a BST-ID H boundary",
        True,
    )


def direct_quantize_points(
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


# ---------------------------------------------------------------------------
# Dense index-space shell/object construction
# ---------------------------------------------------------------------------

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


def close_shell(shell: np.ndarray, iterations: int):
    if iterations <= 0:
        return shell.copy()

    structure = np.ones((3, 3, 3), dtype=bool)
    return ndimage.binary_closing(
        shell,
        structure=structure,
        iterations=iterations,
    )


def directional_fill(shell: np.ndarray, require_ns: bool):
    # x index increases eastward.
    west = np.zeros_like(shell)
    c = np.maximum.accumulate(shell, axis=0)
    west[1:, :, :] = c[:-1, :, :]

    east = np.zeros_like(shell)
    c = np.maximum.accumulate(shell[::-1, :, :], axis=0)[::-1, :, :]
    east[:-1, :, :] = c[1:, :, :]

    # H index increases with altitude.
    up = np.zeros_like(shell)
    c = np.maximum.accumulate(shell[:, :, ::-1], axis=2)[:, :, ::-1]
    up[:, :, :-1] = c[:, :, 1:]

    interior = (~shell) & west & east & up

    if require_ns:
        # Slippy-map y increases southward, but enclosure is symmetric.
        one_side = np.zeros_like(shell)
        c = np.maximum.accumulate(shell, axis=1)
        one_side[:, 1:, :] = c[:, :-1, :]

        other_side = np.zeros_like(shell)
        c = np.maximum.accumulate(shell[:, ::-1, :], axis=1)[:, ::-1, :]
        other_side[:, :-1, :] = c[:, 1:, :]

        interior &= one_side & other_side

    return shell | interior, interior


def flood_fill(shell: np.ndarray):
    structure = ndimage.generate_binary_structure(3, 1)
    solid = ndimage.binary_fill_holes(shell, structure=structure)
    return solid, solid & (~shell)


# ---------------------------------------------------------------------------
# BST-ID cells / normalization
# ---------------------------------------------------------------------------

def triplets_to_cells(triplets: np.ndarray, zx: int, zy: int, zh: int):
    return [
        Cell.from_parts(
            x=(zx, int(ix)),
            y=(zy, int(iy)),
            f=(zh, int(ih)),
        )
        for ix, iy, ih in triplets
    ]


def normalize_cells(cells, mode: str):
    if mode == "indexed-batch":
        if not HAVE_INDEXED_BATCH:
            raise RuntimeError(
                "normalize_indexed_batch() was not found in bst_id.region_index. "
                "Use the indexed implementation bundle or --normalizer reference."
            )
        return normalize_indexed_batch(cells)

    if mode == "reference":
        return normalize(cells)

    raise ValueError(mode)


def cell_to_row(c):
    z = tuple(c.zooms)
    v = tuple(c.values)
    return {
        "zx": int(z[0]), "ix": int(v[0]),
        "zy": int(z[1]), "iy": int(v[1]),
        "zh": int(z[2]), "ih": int(v[2]),
    }


# ---------------------------------------------------------------------------
# Prefix bounds / local centers
# ---------------------------------------------------------------------------

def x_bounds(z: int, i: int):
    n = float(1 << z)
    return i / n * 360.0 - 180.0, (i + 1) / n * 360.0 - 180.0


def y_to_lat(y: float, z: int):
    n = float(1 << z)
    v = math.pi * (1.0 - 2.0 * y / n)
    return math.degrees(math.atan(math.sinh(v)))


def y_bounds(z: int, i: int):
    north = y_to_lat(i, z)
    south = y_to_lat(i + 1, z)
    return south, north


def h_bounds(z: int, i: int):
    n = float(1 << z)
    return (
        H_MIN + (i / n) * H_SPAN,
        H_MIN + ((i + 1) / n) * H_SPAN,
    )


def row_bounds(row: dict, height_origin: float):
    west, east = x_bounds(row["zx"], row["ix"])
    south, north = y_bounds(row["zy"], row["iy"])
    h0, h1 = h_bounds(row["zh"], row["ih"])
    return {
        "west": west,
        "south": south,
        "east": east,
        "north": north,
        "bottom_u": h0 - height_origin,
        "top_u": h1 - height_origin,
    }


def indices_to_local_centers(
    indices: np.ndarray,
    zx: int,
    zy: int,
    zh: int,
    anchor_lat: float,
    anchor_lon: float,
    height_origin: float,
):
    to_local = make_wgs84_to_local(anchor_lat, anchor_lon)

    n_x = float(1 << zx)
    lon = ((indices[:, 0] + 0.5) / n_x) * 360.0 - 180.0

    n_y = float(1 << zy)
    y = indices[:, 1].astype(np.float64) + 0.5
    lat = np.array([y_to_lat(v, zy) for v in y])

    east, north = to_local.transform(lon, lat)

    width_h = H_SPAN / float(1 << zh)
    h = H_MIN + (indices[:, 2] + 0.5) * width_h
    up = h - height_origin

    return np.column_stack((east, north, up))


def save_ply_centers(path: Path, centers: np.ndarray):
    p = o3d.geometry.PointCloud()
    p.points = o3d.utility.Vector3dVector(centers)
    if not o3d.io.write_point_cloud(str(path), p, write_ascii=False):
        raise RuntimeError(f"Could not save {path}")


# ---------------------------------------------------------------------------
# Maximum inscribed prefix
# ---------------------------------------------------------------------------

def integral_volume(occ):
    sat = occ.astype(np.int32).cumsum(0).cumsum(1).cumsum(2)
    return np.pad(sat, ((1, 0), (1, 0), (1, 0)), mode="constant")


def box_sum(sat, x0, x1, y0, y1, z0, z1):
    return (
        sat[x1, y1, z1]
        - sat[x0, y1, z1]
        - sat[x1, y0, z1]
        - sat[x1, y1, z0]
        + sat[x0, y0, z1]
        + sat[x0, y1, z0]
        + sat[x1, y0, z0]
        - sat[x0, y0, z0]
    )


def max_inscribed_prefix(fine: np.ndarray, zx: int, zy: int, zh: int):
    mins = fine.min(axis=0).astype(np.int64)
    maxs = fine.max(axis=0).astype(np.int64)
    shape = (maxs - mins + 1).astype(int)

    occ = np.zeros(tuple(shape), dtype=bool)
    loc = fine - mins
    occ[loc[:, 0], loc[:, 1], loc[:, 2]] = True
    sat = integral_volume(occ)

    mdx = min(zx - 1, int(math.floor(math.log2(shape[0]))))
    mdy = min(zy - 1, int(math.floor(math.log2(shape[1]))))
    mdh = min(zh - 1, int(math.floor(math.log2(shape[2]))))

    def starts(vmin, vmax, block):
        first = ((int(vmin) + block - 1) // block) * block
        last = ((int(vmax) + 1 - block) // block) * block
        if last < first:
            return ()
        return range(first, last + 1, block)

    for s in range(mdx + mdy + mdh, -1, -1):
        found = []
        for dx in range(mdx + 1):
            for dy in range(mdy + 1):
                dh = s - dx - dy
                if dh < 0 or dh > mdh:
                    continue

                bx, by, bh = 1 << dx, 1 << dy, 1 << dh
                cap = bx * by * bh

                for sx in starts(mins[0], maxs[0], bx):
                    lx0, lx1 = sx - mins[0], sx - mins[0] + bx
                    for sy in starts(mins[1], maxs[1], by):
                        ly0, ly1 = sy - mins[1], sy - mins[1] + by
                        for sh in starts(mins[2], maxs[2], bh):
                            lz0, lz1 = sh - mins[2], sh - mins[2] + bh
                            if box_sum(sat, lx0, lx1, ly0, ly1, lz0, lz1) == cap:
                                found.append({
                                    "zx": zx - dx, "ix": sx >> dx,
                                    "zy": zy - dy, "iy": sy >> dy,
                                    "zh": zh - dh, "ih": sh >> dh,
                                    "fine_atom_capacity": cap,
                                    "coarsening_bits": {"x": dx, "y": dy, "h": dh},
                                })
        if found:
            return found

    raise RuntimeError("No maximum inscribed cell found.")


def metric_volume(row: dict, height_origin: float, anchor_lat: float):
    b = row_bounds(row, height_origin)
    lat_mid = 0.5 * (b["south"] + b["north"])
    # Local Earth-scale approximation is adequate only for tie-breaking.
    r = 6378137.0
    dx = math.radians(b["east"] - b["west"]) * r * math.cos(math.radians(lat_mid))
    dy = math.radians(b["north"] - b["south"]) * r
    dz = b["top_u"] - b["bottom_u"]
    return abs(dx * dy * dz)


# ---------------------------------------------------------------------------
# Minimum enclosing cell / cover
# ---------------------------------------------------------------------------

def lcp_axis(vmin: int, vmax: int, fine_z: int):
    x = int(vmin) ^ int(vmax)
    q = fine_z if x == 0 else fine_z - x.bit_length()
    if q < 1:
        return None, None
    return q, int(vmin) >> (fine_z - q)


def enclosing_for_subset(sub: np.ndarray, zx: int, zy: int, zh: int):
    mn, mx = sub.min(axis=0), sub.max(axis=0)
    qx, ix = lcp_axis(mn[0], mx[0], zx)
    qy, iy = lcp_axis(mn[1], mx[1], zy)
    qh, ih = lcp_axis(mn[2], mx[2], zh)
    if qx is None or qy is None or qh is None:
        return None
    return {"zx": qx, "ix": ix, "zy": qy, "iy": iy, "zh": qh, "ih": ih}


def minimum_enclosing(fine: np.ndarray, zx: int, zy: int, zh: int):
    single = enclosing_for_subset(fine, zx, zy, zh)
    if single is not None:
        return "single_cell", [single]

    root_bits = np.column_stack((
        fine[:, 0] >> (zx - 1),
        fine[:, 1] >> (zy - 1),
        fine[:, 2] >> (zh - 1),
    ))

    cover = []
    for combo in np.unique(root_bits, axis=0):
        mask = np.all(root_bits == combo, axis=1)
        c = enclosing_for_subset(fine[mask], zx, zy, zh)
        if c is None:
            raise RuntimeError("Unexpected invalid root-partition cover.")
        c["root_bits"] = {
            "x": int(combo[0]),
            "y": int(combo[1]),
            "h": int(combo[2]),
        }
        cover.append(c)

    return "minimum_valid_cover", cover


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def save_cells_csv(path: Path, cells):
    with path.open("w", newline="", encoding="utf-8") as f:
        fields = ["zx", "ix", "zy", "iy", "zh", "ih"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for c in cells:
            w.writerow(cell_to_row(c) if not isinstance(c, dict) else c)


def zoom_stats(cells):
    ctr = Counter()
    for c in cells:
        r = cell_to_row(c)
        ctr[(r["zx"], r["zy"], r["zh"])] += 1
    return [
        {"zx": z[0], "zy": z[1], "zh": z[2], "count": n}
        for z, n in ctr.most_common()
    ]


def verify_exact(fine_cells, normalized, zx, zy, zh):
    wz = {"x": zx, "y": zy, "f": zh}
    a = region_atoms(fine_cells, wz)
    b = region_atoms(normalized, wz)
    return a == b, len(a), len(b)


def cesium_records(rows, height_origin: float):
    return [
        {
            **r,
            **row_bounds(r, height_origin),
        }
        for r in rows
    ]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("ply", type=Path)
    ap.add_argument(
        "--transform-json",
        type=Path,
        default=Path("tsuruga_enu_anchor_transform.json"),
    )

    ap.add_argument("--min-u", type=float, default=0.0)
    ap.add_argument(
        "--crop",
        type=float,
        nargs=4,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        help="Optional horizontal local EN crop before BST-ID quantization.",
    )

    ap.add_argument("--zx", type=int, default=24)
    ap.add_argument("--zy", type=int, default=24)
    ap.add_argument("--zh", type=int, default=14)

    ap.add_argument("--anchor-alt", type=float, default=None)
    ap.add_argument("--margin-cells", type=int, default=2)

    ap.add_argument(
        "--close-iterations",
        type=int,
        default=1,
    )
    ap.add_argument(
        "--fill-mode",
        choices=["directional", "flood"],
        default="directional",
    )
    ap.add_argument("--require-ns", action="store_true")

    ap.add_argument(
        "--normalizer",
        choices=["indexed-batch", "reference"],
        default="indexed-batch",
    )
    ap.add_argument("--verify", action="store_true")

    ap.add_argument(
        "--out-prefix",
        default="tsuruga_native",
    )

    args = ap.parse_args()
    t_all = time.perf_counter()

    meta = load_transform_metadata(args.transform_json)
    height_origin, height_source, local_height_aligned = resolve_height_origin(
        meta,
        args.zh,
        args.anchor_alt,
    )

    pcd = o3d.io.read_point_cloud(str(args.ply))
    if pcd.is_empty():
        raise RuntimeError(f"Could not load {args.ply}")

    points = np.asarray(pcd.points)
    n_input = len(points)

    keep = points[:, 2] >= args.min_u
    if args.crop is not None:
        xmin, xmax, ymin, ymax = args.crop
        keep &= (
            (points[:, 0] >= xmin) & (points[:, 0] <= xmax) &
            (points[:, 1] >= ymin) & (points[:, 1] <= ymax)
        )
    points = points[keep]

    if len(points) == 0:
        raise RuntimeError("No points remain after crop.")

    print("BST-ID-native shell quantization")
    print("-------------------------------")
    print("input points           :", n_input)
    print("points after crop      :", len(points))
    print("anchor lat/lon         :", meta["lat"], meta["lon"])
    print("fine zoom              :", (args.zx, args.zy, args.zh))
    print("minimum local U [m]    :", args.min_u)
    print("height origin [m]      :", height_origin)
    print("height origin source   :", height_source)

    # 1) Direct point -> fine BST-ID shell.
    t0 = time.perf_counter()
    raw_shell_idx = direct_quantize_points(
        points,
        meta["lat"],
        meta["lon"],
        height_origin,
        args.zx,
        args.zy,
        args.zh,
    )
    print("raw shell BST cells    :", len(raw_shell_idx))
    print("quantization seconds   :", time.perf_counter() - t0)

    # 2) Dense local window in BST integer-index space.
    shell_raw, dense_origin = indices_to_dense(
        raw_shell_idx,
        args.margin_cells,
    )

    shell = close_shell(shell_raw, args.close_iterations)

    # 3) Interior fill IN BST-ID index space.
    if args.fill_mode == "directional":
        solid, interior = directional_fill(shell, args.require_ns)
    else:
        solid, interior = flood_fill(shell)

    shell_idx = dense_to_indices(shell, dense_origin)
    solid_idx = dense_to_indices(solid, dense_origin)

    print()
    print("BST-ID index-space object")
    print("-------------------------")
    print("dense window shape     :", shell.shape)
    print("closed shell cells     :", len(shell_idx))
    print("interior cells         :", int(interior.sum()))
    print("solid fine BST cells   :", len(solid_idx))

    # 4) Fine BST-ID cells -> Normalize.
    fine_cells = triplets_to_cells(solid_idx, args.zx, args.zy, args.zh)

    t0 = time.perf_counter()
    normalized = normalize_cells(fine_cells, args.normalizer)
    norm_s = time.perf_counter() - t0

    n_fine = len(fine_cells)
    n_norm = len(normalized)
    reduction = 1.0 - n_norm / n_fine

    print()
    print("Normalize")
    print("---------")
    print("fine BST-ID cells     :", n_fine)
    print("normalized prefixes   :", n_norm)
    print("reduction             :", f"{100*reduction:.3f}%")
    print("seconds               :", norm_s)

    # 5) Exact verification.
    verification = {
        "requested": bool(args.verify),
        "passed": None,
        "fine_atoms": None,
        "normalized_atoms": None,
    }
    if args.verify:
        ok, na, nb = verify_exact(
            fine_cells, normalized, args.zx, args.zy, args.zh
        )
        verification.update(
            passed=bool(ok),
            fine_atoms=int(na),
            normalized_atoms=int(nb),
        )
        print("exact verification     :", ok)
        if not ok:
            raise RuntimeError("Normalized occupancy mismatch.")

    # 6) Inner/outer special cells.
    inscribed_candidates = max_inscribed_prefix(
        solid_idx, args.zx, args.zy, args.zh
    )
    inscribed = max(
        inscribed_candidates,
        key=lambda r: metric_volume(
            r,
            height_origin,
            meta["lat"],
        ),
    )

    enclosing_mode, enclosing = minimum_enclosing(
        solid_idx,
        args.zx,
        args.zy,
        args.zh,
    )

    # 7) Save local PLY centers for Open3D checks.
    shell_centers = indices_to_local_centers(
        shell_idx,
        args.zx,
        args.zy,
        args.zh,
        meta["lat"],
        meta["lon"],
        height_origin,
    )
    solid_centers = indices_to_local_centers(
        solid_idx,
        args.zx,
        args.zy,
        args.zh,
        meta["lat"],
        meta["lon"],
        height_origin,
    )

    prefix = Path(args.out_prefix)
    shell_ply = Path(f"{prefix}_shell.ply")
    solid_ply = Path(f"{prefix}_solid.ply")
    norm_csv = Path(f"{prefix}_normalized.csv")
    stats_json = Path(f"{prefix}_stats.json")
    cesium_json = Path(f"{prefix}_cesium.json")

    save_ply_centers(shell_ply, shell_centers)
    save_ply_centers(solid_ply, solid_centers)
    save_cells_csv(norm_csv, normalized)

    zstats = zoom_stats(normalized)

    stats = {
        "pipeline": "point cloud -> direct fine BST-ID shell -> BST-ID-index fill -> Normalize",
        "input": {
            "ply": str(args.ply),
            "transform_json": str(args.transform_json),
            "input_points": n_input,
            "points_after_crop": len(points),
            "min_u_m": args.min_u,
            "horizontal_crop": args.crop,
        },
        "georeference": {
            "anchor_latitude_deg": meta["lat"],
            "anchor_longitude_deg": meta["lon"],
            "height_origin_m": height_origin,
            "height_origin_source": height_source,
            "local_u0_aligned_to_h_boundary": local_height_aligned,
        },
        "fine_grid": {
            "zx": args.zx,
            "zy": args.zy,
            "zh": args.zh,
            "raw_shell_cells": len(raw_shell_idx),
            "closed_shell_cells": len(shell_idx),
            "fill_mode": args.fill_mode,
            "require_ns": args.require_ns,
            "close_iterations": args.close_iterations,
            "interior_cells": int(interior.sum()),
            "solid_fine_cells": n_fine,
            "dense_window_shape": list(shell.shape),
        },
        "normalization": {
            "normalizer": args.normalizer,
            "normalized_prefixes": n_norm,
            "reduction_fraction": reduction,
            "reduction_percent": 100.0 * reduction,
            "seconds": norm_s,
            "zoom_vectors": zstats,
        },
        "maximum_inscribed_cell": {
            **inscribed,
            **row_bounds(inscribed, height_origin),
            "metric_volume_m3_approx": metric_volume(
                inscribed,
                height_origin,
                meta["lat"],
            ),
        },
        "minimum_enclosing": {
            "mode": enclosing_mode,
            "cell_count": len(enclosing),
            "cells": [
                {**r, **row_bounds(r, height_origin)}
                for r in enclosing
            ],
        },
        "verification": verification,
        "total_seconds": time.perf_counter() - t_all,
    }

    stats_json.write_text(
        json.dumps(stats, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    # 8) Compact Cesium data.  Keep prefix descriptors plus precomputed bounds.
    fine_rows = [
        {
            "zx": args.zx, "ix": int(v[0]),
            "zy": args.zy, "iy": int(v[1]),
            "zh": args.zh, "ih": int(v[2]),
        }
        for v in solid_idx
    ]
    norm_rows = [cell_to_row(c) for c in normalized]

    cesium = {
        "metadata": {
            "title": "Tsuruga Castle BST-ID object",
            "anchor_latitude_deg": meta["lat"],
            "anchor_longitude_deg": meta["lon"],
            "height_origin_m": height_origin,
            "height_display": "terrain-relative local U",
            "fine_zoom": {"x": args.zx, "y": args.zy, "h": args.zh},
            "fine_cells": n_fine,
            "normalized_prefixes": n_norm,
            "reduction_percent": 100.0 * reduction,
            "source_note": (
                "Tsuruga Castle 3D model; BST-ID-native occupancy generated "
                "from the aligned point cloud."
            ),
        },
        "bitmap": cesium_records(fine_rows, height_origin),
        "normalized": cesium_records(norm_rows, height_origin),
        "maximum_inscribed": {
            **inscribed,
            **row_bounds(inscribed, height_origin),
        },
        "minimum_enclosing": {
            "mode": enclosing_mode,
            "cells": cesium_records(enclosing, height_origin),
        },
    }

    cesium_json.write_text(
        json.dumps(cesium, separators=(",", ":"), ensure_ascii=False),
        encoding="utf-8",
    )

    print()
    print("Special cells")
    print("-------------")
    print(
        "maximum inscribed     :",
        (inscribed["zx"], inscribed["ix"],
         inscribed["zy"], inscribed["iy"],
         inscribed["zh"], inscribed["ih"]),
    )
    print("minimum enclosing mode:", enclosing_mode)
    print("enclosing cell count   :", len(enclosing))

    print()
    print("Saved")
    print("-----")
    print(shell_ply)
    print(solid_ply)
    print(norm_csv)
    print(stats_json)
    print(cesium_json)


if __name__ == "__main__":
    main()
