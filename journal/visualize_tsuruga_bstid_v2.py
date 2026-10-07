#!/usr/bin/env python3
"""
visualize_tsuruga_bstid.py

Visualize, in one Open3D window:

  B : fixed-grid bitmap/voxel occupancy
  H : normalized hierarchical BST-ID fill
  I : toggle maximum inscribed BST-ID cell
  E : toggle minimum enclosing BST-ID cell/cover
  A : show bitmap + hierarchy + both special cells
  C : hide bitmap/hierarchy, keep reference frame/special-cell state
  P : save screenshot

The script also computes:

1. Maximum inscribed BST-ID cell
   The prefix cell with the maximum number of declared fine BST-ID atoms that
   is completely occupied by the object.  This is computed directly from the
   fine occupancy, not merely selected from the normalized output.

2. Minimum enclosing BST-ID cell, or minimum valid cover if one cell is impossible
   The smallest single anisotropic BST-ID prefix cell containing all fine
   occupied cells.  This is obtained from the per-axis longest common prefix
   (LCP) of the fine occupied indices.

Inputs are the outputs of tsuruga_ply_to_bstid_object.py.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
import open3d as o3d
from pyproj import CRS, Transformer


WEB_MERCATOR_LAT_MAX = 85.0511287798066
H_MIN = -16383.0
H_MAX = 16384.0
H_SPAN = H_MAX - H_MIN


# ---------------------------------------------------------------------------
# Coordinate / BST-ID utilities
# ---------------------------------------------------------------------------

def make_local_crs(anchor_lat: float, anchor_lon: float) -> CRS:
    return CRS.from_proj4(
        f"+proj=aeqd +lat_0={anchor_lat:.12f} "
        f"+lon_0={anchor_lon:.12f} +datum=WGS84 +units=m +no_defs"
    )


def local_to_wgs84_transformer(anchor_lat: float, anchor_lon: float):
    return Transformer.from_crs(
        make_local_crs(anchor_lat, anchor_lon),
        CRS.from_epsg(4326),
        always_xy=True,
    )


def wgs84_to_local_transformer(anchor_lat: float, anchor_lon: float):
    return Transformer.from_crs(
        CRS.from_epsg(4326),
        make_local_crs(anchor_lat, anchor_lon),
        always_xy=True,
    )


def lonlat_to_prefix(lon_deg, lat_deg, zx: int, zy: int):
    lon_deg = np.asarray(lon_deg, dtype=np.float64)
    lat_deg = np.asarray(lat_deg, dtype=np.float64)

    lon = ((lon_deg + 180.0) % 360.0) - 180.0
    lat = np.clip(lat_deg, -WEB_MERCATOR_LAT_MAX, WEB_MERCATOR_LAT_MAX)

    nx = float(1 << zx)
    x_norm = (lon + 180.0) / 360.0
    ix = np.floor(x_norm * nx).astype(np.int64)
    ix = np.clip(ix, 0, (1 << zx) - 1)

    ny = float(1 << zy)
    lat_rad = np.radians(lat)
    y_norm = (1.0 - np.arcsinh(np.tan(lat_rad)) / math.pi) / 2.0
    iy = np.floor(y_norm * ny).astype(np.int64)
    iy = np.clip(iy, 0, (1 << zy) - 1)

    return ix, iy


def altitude_to_prefix(h_m, zh: int):
    h_m = np.asarray(h_m, dtype=np.float64)
    n = float(1 << zh)
    u = (h_m - H_MIN) / H_SPAN
    u = np.clip(u, 0.0, np.nextafter(1.0, 0.0))
    return np.floor(u * n).astype(np.int64)


def webmercator_x_bounds(z: int, i: int):
    n = float(1 << z)
    lon0 = i / n * 360.0 - 180.0
    lon1 = (i + 1) / n * 360.0 - 180.0
    return lon0, lon1


def webmercator_y_to_lat(y: float, z: int):
    n = float(1 << z)
    v = math.pi * (1.0 - 2.0 * y / n)
    return math.degrees(math.atan(math.sinh(v)))


def webmercator_y_bounds(z: int, i: int):
    # Slippy-map Y increases southward.
    lat_north = webmercator_y_to_lat(i, z)
    lat_south = webmercator_y_to_lat(i + 1, z)
    return lat_south, lat_north


def altitude_bounds(z: int, i: int):
    n = float(1 << z)
    h0 = H_MIN + (i / n) * H_SPAN
    h1 = H_MIN + ((i + 1) / n) * H_SPAN
    return h0, h1


def prefix_to_local_bounds(
    cell: dict,
    anchor_lat: float,
    anchor_lon: float,
    anchor_alt: float,
    to_local,
):
    zx, ix = int(cell["zx"]), int(cell["ix"])
    zy, iy = int(cell["zy"]), int(cell["iy"])
    zh, ih = int(cell["zh"]), int(cell["ih"])

    lon0, lon1 = webmercator_x_bounds(zx, ix)
    lat0, lat1 = webmercator_y_bounds(zy, iy)

    lons = np.array([lon0, lon0, lon1, lon1])
    lats = np.array([lat0, lat1, lat0, lat1])
    east, north = to_local.transform(lons, lats)

    e0, e1 = float(np.min(east)), float(np.max(east))
    n0, n1 = float(np.min(north)), float(np.max(north))

    h0, h1 = altitude_bounds(zh, ih)
    u0 = h0 - anchor_alt
    u1 = h1 - anchor_alt

    return np.array([e0, n0, u0]), np.array([e1, n1, u1])


def encode_fine_cells(
    centers_enu: np.ndarray,
    anchor_lat: float,
    anchor_lon: float,
    anchor_alt: float,
    zx: int,
    zy: int,
    zh: int,
):
    to_wgs84 = local_to_wgs84_transformer(anchor_lat, anchor_lon)

    lon, lat = to_wgs84.transform(
        centers_enu[:, 0],
        centers_enu[:, 1],
    )

    ix, iy = lonlat_to_prefix(lon, lat, zx, zy)
    ih = altitude_to_prefix(anchor_alt + centers_enu[:, 2], zh)

    triplets = np.column_stack((ix, iy, ih))
    return np.unique(triplets, axis=0)


# ---------------------------------------------------------------------------
# Exact discrete maximum-inscribed cell
# ---------------------------------------------------------------------------

def integral_volume(occ: np.ndarray):
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
    """
    Find the occupied BST-ID cell containing the largest number of fine atoms.

    Search is performed on globally aligned dyadic blocks using a summed-volume
    table.  Among candidates with the same maximum fine-atom capacity, the
    caller can tie-break by physical metric volume.
    """
    mins = fine.min(axis=0).astype(np.int64)
    maxs = fine.max(axis=0).astype(np.int64)
    shape = (maxs - mins + 1).astype(int)

    occ = np.zeros(tuple(shape), dtype=bool)
    local = fine - mins
    occ[local[:, 0], local[:, 1], local[:, 2]] = True
    sat = integral_volume(occ)

    max_dx = min(zx - 1, int(math.floor(math.log2(shape[0]))))
    max_dy = min(zy - 1, int(math.floor(math.log2(shape[1]))))
    max_dh = min(zh - 1, int(math.floor(math.log2(shape[2]))))

    def aligned_starts(vmin, vmax, block):
        first = ((int(vmin) + block - 1) // block) * block
        last = ((int(vmax) + 1 - block) // block) * block
        if last < first:
            return ()
        return range(first, last + 1, block)

    max_s = max_dx + max_dy + max_dh

    for s in range(max_s, -1, -1):
        found = []

        for dx in range(max_dx + 1):
            for dy in range(max_dy + 1):
                dh = s - dx - dy
                if dh < 0 or dh > max_dh:
                    continue

                bx = 1 << dx
                by = 1 << dy
                bh = 1 << dh
                capacity = bx * by * bh

                xs = aligned_starts(mins[0], maxs[0], bx)
                ys = aligned_starts(mins[1], maxs[1], by)
                zs = aligned_starts(mins[2], maxs[2], bh)

                for sx in xs:
                    lx0 = sx - mins[0]
                    lx1 = lx0 + bx
                    for sy in ys:
                        ly0 = sy - mins[1]
                        ly1 = ly0 + by
                        for sh in zs:
                            lz0 = sh - mins[2]
                            lz1 = lz0 + bh

                            if box_sum(
                                sat,
                                lx0, lx1,
                                ly0, ly1,
                                lz0, lz1,
                            ) == capacity:
                                found.append(
                                    {
                                        "zx": zx - dx,
                                        "ix": sx >> dx,
                                        "zy": zy - dy,
                                        "iy": sy >> dy,
                                        "zh": zh - dh,
                                        "ih": sh >> dh,
                                        "fine_atom_capacity": capacity,
                                        "coarsening_bits": {
                                            "x": dx,
                                            "y": dy,
                                            "h": dh,
                                        },
                                    }
                                )

        if found:
            return found

    raise RuntimeError("No inscribed BST-ID cell found.")


# ---------------------------------------------------------------------------
# Minimum enclosing cell
# ---------------------------------------------------------------------------

def lcp_cell(vmin: int, vmax: int, fine_z: int, allow_zero: bool = False):
    """Return the longest common-prefix cell on one axis."""
    x = int(vmin) ^ int(vmax)
    common = fine_z if x == 0 else fine_z - x.bit_length()

    if common < 1 and not allow_zero:
        return None, None

    if common == 0:
        return 0, 0

    prefix = int(vmin) >> (fine_z - common)
    return common, prefix


def _enclosing_cell_for_subset(
    subset: np.ndarray,
    zx: int,
    zy: int,
    zh: int,
):
    """Smallest single valid BST-ID cell containing one subset."""
    mins = subset.min(axis=0)
    maxs = subset.max(axis=0)

    qx, ix = lcp_cell(mins[0], maxs[0], zx)
    qy, iy = lcp_cell(mins[1], maxs[1], zy)
    qh, ih = lcp_cell(mins[2], maxs[2], zh)

    if qx is None or qy is None or qh is None:
        return None

    return {
        "zx": qx,
        "ix": ix,
        "zy": qy,
        "iy": iy,
        "zh": qh,
        "ih": ih,
        "fine_atom_capacity": (
            (1 << (zx - qx))
            * (1 << (zy - qy))
            * (1 << (zh - qh))
        ),
        "coarsening_bits": {
            "x": zx - qx,
            "y": zy - qy,
            "h": zh - qh,
        },
    }


def minimum_enclosing_prefix(fine: np.ndarray, zx: int, zy: int, zh: int):
    """
    Smallest SINGLE valid BST-ID cell containing all occupied fine cells.

    Returns None if the object crosses a zoom-1 root split on any active axis.
    Since the current BST-ID definition has minimum zoom 1, such a cell does
    not exist inside the valid BST-ID domain.
    """
    return _enclosing_cell_for_subset(fine, zx, zy, zh)


def minimum_enclosing_cover(fine: np.ndarray, zx: int, zy: int, zh: int):
    """
    Minimum valid enclosing COVER when no single BST-ID cell exists.

    A valid BST-ID cell cannot cross a zoom-1 split on any active axis.
    Therefore the occupied fine cells are partitioned by their first
    (zoom-1) bit on x/y/h.  Each occupied root-bit combination requires at
    least one valid cell, and the LCP cell of that subset is the smallest
    single cell enclosing that subset.

    The returned cover is therefore minimal in the number of zoom>=1 root
    partitions and tight within each occupied partition.
    """
    root_bits = np.column_stack(
        (
            fine[:, 0] >> (zx - 1),
            fine[:, 1] >> (zy - 1),
            fine[:, 2] >> (zh - 1),
        )
    )

    combos = np.unique(root_bits, axis=0)
    cover = []

    for combo in combos:
        mask = np.all(root_bits == combo, axis=1)
        subset = fine[mask]

        cell = _enclosing_cell_for_subset(subset, zx, zy, zh)
        if cell is None:
            raise RuntimeError(
                "Internal error: a root-partition subset still has no "
                "valid zoom>=1 enclosing cell."
            )

        cell["root_bits"] = {
            "x": int(combo[0]),
            "y": int(combo[1]),
            "h": int(combo[2]),
        }
        cell["occupied_fine_cells_in_partition"] = int(len(subset))
        cover.append(cell)

    return cover


# ---------------------------------------------------------------------------
# Geometry creation
# ---------------------------------------------------------------------------

BOX_TRIANGLES = np.array(
    [
        [0, 1, 2], [0, 2, 3],  # bottom
        [4, 6, 5], [4, 7, 6],  # top
        [0, 4, 5], [0, 5, 1],
        [1, 5, 6], [1, 6, 2],
        [2, 6, 7], [2, 7, 3],
        [3, 7, 4], [3, 4, 0],
    ],
    dtype=np.int32,
)

BOX_EDGES = np.array(
    [
        [0, 1], [1, 2], [2, 3], [3, 0],
        [4, 5], [5, 6], [6, 7], [7, 4],
        [0, 4], [1, 5], [2, 6], [3, 7],
    ],
    dtype=np.int32,
)


def box_vertices(bmin, bmax, shrink=1.0):
    bmin = np.asarray(bmin, dtype=float)
    bmax = np.asarray(bmax, dtype=float)

    if shrink < 1.0:
        c = 0.5 * (bmin + bmax)
        half = 0.5 * (bmax - bmin) * shrink
        bmin = c - half
        bmax = c + half

    x0, y0, z0 = bmin
    x1, y1, z1 = bmax

    return np.array(
        [
            [x0, y0, z0],
            [x1, y0, z0],
            [x1, y1, z0],
            [x0, y1, z0],
            [x0, y0, z1],
            [x1, y0, z1],
            [x1, y1, z1],
            [x0, y1, z1],
        ],
        dtype=float,
    )


def boxes_to_mesh(bounds, color, shrink=1.0):
    vertices = []
    triangles = []

    for bmin, bmax in bounds:
        base = len(vertices)
        vertices.extend(box_vertices(bmin, bmax, shrink=shrink))
        triangles.extend(BOX_TRIANGLES + base)

    mesh = o3d.geometry.TriangleMesh(
        o3d.utility.Vector3dVector(np.asarray(vertices)),
        o3d.utility.Vector3iVector(np.asarray(triangles)),
    )
    mesh.compute_vertex_normals()
    mesh.paint_uniform_color(color)
    return mesh


def boxes_to_lines(bounds, color, shrink=1.0):
    points = []
    lines = []

    for bmin, bmax in bounds:
        base = len(points)
        points.extend(box_vertices(bmin, bmax, shrink=shrink))
        lines.extend(BOX_EDGES + base)

    ls = o3d.geometry.LineSet(
        points=o3d.utility.Vector3dVector(np.asarray(points)),
        lines=o3d.utility.Vector2iVector(np.asarray(lines)),
    )
    ls.colors = o3d.utility.Vector3dVector(
        np.tile(np.asarray(color, dtype=float), (len(lines), 1))
    )
    return ls


def make_anchor_marker(radius=1.0):
    sphere = o3d.geometry.TriangleMesh.create_sphere(radius=radius)
    sphere.paint_uniform_color([0.95, 0.75, 0.1])
    sphere.translate([0.0, 0.0, 0.0])
    sphere.compute_vertex_normals()
    return sphere


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_normalized_csv(path: Path):
    cells = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cells.append(
                {
                    "zx": int(row["zx"]),
                    "ix": int(row["ix"]),
                    "zy": int(row["zy"]),
                    "iy": int(row["iy"]),
                    "zh": int(row["zh"]),
                    "ih": int(row["ih"]),
                }
            )
    return cells


def metric_volume(bounds):
    bmin, bmax = bounds
    d = np.maximum(np.asarray(bmax) - np.asarray(bmin), 0.0)
    return float(d[0] * d[1] * d[2])


# ---------------------------------------------------------------------------
# Visualizer
# ---------------------------------------------------------------------------

class ToggleViewer:
    def __init__(
        self,
        bitmap,
        hierarchy_mesh,
        hierarchy_lines,
        inscribed_mesh,
        inscribed_lines,
        enclosing_lines,
        frame,
        anchor_marker,
        screenshot_path,
    ):
        self.bitmap = bitmap
        self.hierarchy_mesh = hierarchy_mesh
        self.hierarchy_lines = hierarchy_lines
        self.inscribed_mesh = inscribed_mesh
        self.inscribed_lines = inscribed_lines
        self.enclosing_lines = enclosing_lines
        self.frame = frame
        self.anchor_marker = anchor_marker
        self.screenshot_path = screenshot_path

        self.vis = o3d.visualization.VisualizerWithKeyCallback()
        self.visible = set()

    def _add(self, name, geom, reset=False):
        if name in self.visible:
            return
        try:
            self.vis.add_geometry(geom, reset_bounding_box=reset)
        except TypeError:
            self.vis.add_geometry(geom)
        self.visible.add(name)

    def _remove(self, name, geom):
        if name not in self.visible:
            return
        try:
            self.vis.remove_geometry(geom, reset_bounding_box=False)
        except TypeError:
            self.vis.remove_geometry(geom)
        self.visible.discard(name)

    def bitmap_only(self, vis):
        self._remove("hier_mesh", self.hierarchy_mesh)
        self._remove("hier_lines", self.hierarchy_lines)
        self._add("bitmap", self.bitmap)
        return False

    def hierarchy_only(self, vis):
        self._remove("bitmap", self.bitmap)
        self._add("hier_mesh", self.hierarchy_mesh)
        self._add("hier_lines", self.hierarchy_lines)
        return False

    def toggle_inscribed(self, vis):
        if "inscribed_mesh" in self.visible:
            self._remove("inscribed_mesh", self.inscribed_mesh)
            self._remove("inscribed_lines", self.inscribed_lines)
        else:
            self._add("inscribed_mesh", self.inscribed_mesh)
            self._add("inscribed_lines", self.inscribed_lines)
        return False

    def toggle_enclosing(self, vis):
        if "enclosing_lines" in self.visible:
            self._remove("enclosing_lines", self.enclosing_lines)
        else:
            self._add("enclosing_lines", self.enclosing_lines)
        return False

    def show_all(self, vis):
        self._add("bitmap", self.bitmap)
        self._add("hier_mesh", self.hierarchy_mesh)
        self._add("hier_lines", self.hierarchy_lines)
        self._add("inscribed_mesh", self.inscribed_mesh)
        self._add("inscribed_lines", self.inscribed_lines)
        self._add("enclosing_lines", self.enclosing_lines)
        return False

    def clear_representation(self, vis):
        self._remove("bitmap", self.bitmap)
        self._remove("hier_mesh", self.hierarchy_mesh)
        self._remove("hier_lines", self.hierarchy_lines)
        return False

    def screenshot(self, vis):
        vis.capture_screen_image(str(self.screenshot_path), do_render=True)
        print("Saved screenshot:", self.screenshot_path)
        return False

    def run(self):
        self.vis.create_window(
            window_name="Tsuruga Castle: Bitmap vs normalized BST-ID",
            width=1500,
            height=950,
        )

        # Permanent reference geometries.
        self._add("frame", self.frame, reset=True)
        self._add("anchor", self.anchor_marker)

        # Initial view: hierarchical representation + both special cells.
        self._add("hier_mesh", self.hierarchy_mesh)
        self._add("hier_lines", self.hierarchy_lines)
        self._add("inscribed_mesh", self.inscribed_mesh)
        self._add("inscribed_lines", self.inscribed_lines)
        self._add("enclosing_lines", self.enclosing_lines)

        self.vis.register_key_callback(ord("B"), self.bitmap_only)
        self.vis.register_key_callback(ord("H"), self.hierarchy_only)
        self.vis.register_key_callback(ord("I"), self.toggle_inscribed)
        self.vis.register_key_callback(ord("E"), self.toggle_enclosing)
        self.vis.register_key_callback(ord("A"), self.show_all)
        self.vis.register_key_callback(ord("C"), self.clear_representation)
        self.vis.register_key_callback(ord("P"), self.screenshot)

        print()
        print("Viewer keys")
        print("-----------")
        print("B : bitmap/fixed-grid occupancy")
        print("H : normalized hierarchical BST-ID fill")
        print("I : toggle maximum inscribed BST-ID cell")
        print("E : toggle minimum enclosing BST-ID cell")
        print("A : show all")
        print("C : hide bitmap/hierarchy")
        print("P : save screenshot")
        print("Q / Esc : close")

        self.vis.run()
        self.vis.destroy_window()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--solid-ply",
        type=Path,
        default=Path("tsuruga_1m_solid_voxels.ply"),
    )
    ap.add_argument(
        "--normalized-csv",
        type=Path,
        default=Path("tsuruga_1m_normalized_prefixes.csv"),
    )
    ap.add_argument(
        "--stats-json",
        type=Path,
        default=Path("tsuruga_1m_stats.json"),
    )

    ap.add_argument(
        "--hier-shrink",
        type=float,
        default=0.96,
        help="Visual shrink factor for normalized boxes to expose cell boundaries",
    )
    ap.add_argument(
        "--bitmap-color",
        nargs=3,
        type=float,
        default=(0.72, 0.72, 0.72),
    )
    ap.add_argument(
        "--hier-color",
        nargs=3,
        type=float,
        default=(0.25, 0.55, 0.90),
    )
    ap.add_argument(
        "--inscribed-color",
        nargs=3,
        type=float,
        default=(0.15, 0.85, 0.25),
    )
    ap.add_argument(
        "--enclosing-color",
        nargs=3,
        type=float,
        default=(0.95, 0.15, 0.15),
    )

    ap.add_argument(
        "--out-report",
        type=Path,
        default=Path("tsuruga_bstid_visualizer_report.json"),
    )
    ap.add_argument(
        "--screenshot",
        type=Path,
        default=Path("tsuruga_bstid_visualizer.png"),
    )
    ap.add_argument(
        "--export-hierarchy-obj",
        type=Path,
        default=Path("tsuruga_hierarchical_bstid.obj"),
    )

    args = ap.parse_args()

    stats = json.loads(args.stats_json.read_text(encoding="utf-8"))
    fine_zoom = stats["bstid"]["fine_zoom"]
    zx = int(fine_zoom["x"])
    zy = int(fine_zoom["y"])
    zh = int(fine_zoom["h"])

    geo = stats["georeference"]
    anchor_lat = float(geo["anchor_latitude_deg"])
    anchor_lon = float(geo["anchor_longitude_deg"])
    anchor_alt = float(geo["anchor_altitude_m"])
    voxel_size = float(stats["local_grid"]["voxel_size_m"])

    normalized = read_normalized_csv(args.normalized_csv)

    pcd = o3d.io.read_point_cloud(str(args.solid_ply))
    if pcd.is_empty():
        raise RuntimeError(f"Could not load {args.solid_ply}")

    centers = np.asarray(pcd.points)

    # ------------------------------------------------------------------
    # Fine BST-ID occupancy
    # ------------------------------------------------------------------
    fine = encode_fine_cells(
        centers,
        anchor_lat,
        anchor_lon,
        anchor_alt,
        zx,
        zy,
        zh,
    )

    expected_fine = int(stats["bstid"]["unique_fine_cells"])
    if len(fine) != expected_fine:
        print(
            "WARNING: re-encoded fine BST-ID count differs from stats:",
            len(fine),
            "vs",
            expected_fine,
        )

    # ------------------------------------------------------------------
    # Special cells
    # ------------------------------------------------------------------
    to_local = wgs84_to_local_transformer(anchor_lat, anchor_lon)

    inscribed_candidates = max_inscribed_prefix(fine, zx, zy, zh)

    best_inscribed = None
    best_inscribed_bounds = None
    best_inscribed_volume = -1.0

    for c in inscribed_candidates:
        bounds = prefix_to_local_bounds(
            c,
            anchor_lat,
            anchor_lon,
            anchor_alt,
            to_local,
        )
        vol = metric_volume(bounds)
        if vol > best_inscribed_volume:
            best_inscribed = c
            best_inscribed_bounds = bounds
            best_inscribed_volume = vol

    enclosing_single = minimum_enclosing_prefix(fine, zx, zy, zh)

    if enclosing_single is not None:
        enclosing_cover = [enclosing_single]
        enclosing_mode = "single_cell"
    else:
        enclosing_cover = minimum_enclosing_cover(fine, zx, zy, zh)
        enclosing_mode = "minimum_valid_cover"
        print()
        print(
            "NOTE: no single valid zoom>=1 BST-ID cell contains the object."
        )
        print(
            "      The object crosses at least one zoom-1 root split."
        )
        print(
            "      Visualizing the minimum valid enclosing BST-ID cover instead:",
            len(enclosing_cover),
            "cells",
        )

    enclosing_bounds_list = [
        prefix_to_local_bounds(
            c,
            anchor_lat,
            anchor_lon,
            anchor_alt,
            to_local,
        )
        for c in enclosing_cover
    ]
    enclosing_volume = sum(metric_volume(b) for b in enclosing_bounds_list)
    enclosing_capacity = sum(
        int(c["fine_atom_capacity"]) for c in enclosing_cover
    )

    # ------------------------------------------------------------------
    # Normalized hierarchy bounds
    # ------------------------------------------------------------------
    hierarchy_bounds = [
        prefix_to_local_bounds(
            c,
            anchor_lat,
            anchor_lon,
            anchor_alt,
            to_local,
        )
        for c in normalized
    ]

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------
    pcd.paint_uniform_color(args.bitmap_color)
    bitmap = o3d.geometry.VoxelGrid.create_from_point_cloud(
        pcd,
        voxel_size=voxel_size,
    )

    hierarchy_mesh = boxes_to_mesh(
        hierarchy_bounds,
        color=args.hier_color,
        shrink=args.hier_shrink,
    )
    hierarchy_lines = boxes_to_lines(
        hierarchy_bounds,
        color=[0.10, 0.10, 0.10],
        shrink=args.hier_shrink,
    )

    inscribed_mesh = boxes_to_mesh(
        [best_inscribed_bounds],
        color=args.inscribed_color,
        shrink=0.98,
    )
    inscribed_lines = boxes_to_lines(
        [best_inscribed_bounds],
        color=[0.0, 0.25, 0.0],
        shrink=0.98,
    )

    enclosing_lines = boxes_to_lines(
        enclosing_bounds_list,
        color=args.enclosing_color,
        shrink=1.0,
    )

    frame_size = max(10.0, 0.08 * np.max(pcd.get_axis_aligned_bounding_box().get_extent()))
    frame = o3d.geometry.TriangleMesh.create_coordinate_frame(
        size=frame_size,
        origin=[0.0, 0.0, 0.0],
    )
    anchor_marker = make_anchor_marker(radius=max(0.5, voxel_size * 0.6))

    # Windows 3D viewers often handle triangle OBJ more reliably than
    # vertex-only point-cloud PLY.  Export the hierarchy as a mesh OBJ.
    o3d.io.write_triangle_mesh(
        str(args.export_hierarchy_obj),
        hierarchy_mesh,
        write_triangle_uvs=False,
    )

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------
    object_fine_cells = len(fine)
    inscribed_capacity = int(best_inscribed["fine_atom_capacity"])
    report = {
        "representation": {
            "bitmap_fine_cells": object_fine_cells,
            "normalized_prefixes": len(normalized),
            "normalization_reduction_percent":
                100.0 * (1.0 - len(normalized) / object_fine_cells),
        },
        "fine_zoom": {
            "x": zx,
            "y": zy,
            "h": zh,
        },
        "maximum_inscribed_bstid_cell": {
            **best_inscribed,
            "local_min_m": best_inscribed_bounds[0].tolist(),
            "local_max_m": best_inscribed_bounds[1].tolist(),
            "metric_volume_m3": best_inscribed_volume,
            "fraction_of_object_fine_atoms":
                inscribed_capacity / object_fine_cells,
        },
        "minimum_enclosing_bstid": {
            "mode": enclosing_mode,
            "single_cell_available": enclosing_single is not None,
            "cell_count": len(enclosing_cover),
            "cells": [
                {
                    **cell,
                    "local_min_m": bounds[0].tolist(),
                    "local_max_m": bounds[1].tolist(),
                    "metric_volume_m3": metric_volume(bounds),
                }
                for cell, bounds in zip(
                    enclosing_cover,
                    enclosing_bounds_list,
                )
            ],
            "total_fine_atom_capacity": int(enclosing_capacity),
            "total_metric_volume_m3": float(enclosing_volume),
            "object_occupancy_fraction_of_enclosing_fine_capacity":
                object_fine_cells / enclosing_capacity,
        },
        "visualization": {
            "bitmap": "Open3D VoxelGrid from solid fine-grid centers",
            "hierarchy": "Filled boxes from normalized BST-ID prefixes",
            "hierarchy_visual_shrink": args.hier_shrink,
            "keys": {
                "B": "bitmap only",
                "H": "hierarchical only",
                "I": "toggle maximum inscribed cell",
                "E": "toggle minimum enclosing cell/cover",
                "A": "show all",
                "C": "hide bitmap/hierarchy",
                "P": "save screenshot",
            },
        },
        "outputs": {
            "report_json": str(args.out_report),
            "hierarchy_obj": str(args.export_hierarchy_obj),
            "screenshot_path": str(args.screenshot),
        },
    }

    args.out_report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("Representation")
    print("--------------")
    print("bitmap fine cells     :", object_fine_cells)
    print("normalized prefixes   :", len(normalized))
    print(
        "reduction             :",
        f"{report['representation']['normalization_reduction_percent']:.3f}%",
    )

    print()
    print("Maximum inscribed BST-ID cell")
    print("------------------------------")
    print(
        "prefix:",
        (
            best_inscribed["zx"], best_inscribed["ix"],
            best_inscribed["zy"], best_inscribed["iy"],
            best_inscribed["zh"], best_inscribed["ih"],
        ),
    )
    print("fine-atom capacity    :", inscribed_capacity)
    print("metric volume [m^3]   :", f"{best_inscribed_volume:.3f}")
    print("local min [m]         :", best_inscribed_bounds[0])
    print("local max [m]         :", best_inscribed_bounds[1])

    print()
    if enclosing_single is not None:
        print("Minimum enclosing BST-ID cell")
        print("------------------------------")
    else:
        print("Minimum valid enclosing BST-ID cover")
        print("------------------------------------")

    print("mode                  :", enclosing_mode)
    print("number of cells       :", len(enclosing_cover))
    print("total fine capacity   :", enclosing_capacity)
    print("total metric volume   :", f"{enclosing_volume:.3f} m^3")

    for idx, (cell, bounds) in enumerate(
        zip(enclosing_cover, enclosing_bounds_list),
        start=1,
    ):
        print(
            f"  [{idx}] prefix:",
            (
                cell["zx"], cell["ix"],
                cell["zy"], cell["iy"],
                cell["zh"], cell["ih"],
            ),
        )
        if "root_bits" in cell:
            print("      root bits       :", cell["root_bits"])
        print("      local min [m]   :", bounds[0])
        print("      local max [m]   :", bounds[1])

    print()
    print("Saved report          :", args.out_report)
    print("Saved hierarchy OBJ   :", args.export_hierarchy_obj)

    # ------------------------------------------------------------------
    # Interactive viewer
    # ------------------------------------------------------------------
    viewer = ToggleViewer(
        bitmap=bitmap,
        hierarchy_mesh=hierarchy_mesh,
        hierarchy_lines=hierarchy_lines,
        inscribed_mesh=inscribed_mesh,
        inscribed_lines=inscribed_lines,
        enclosing_lines=enclosing_lines,
        frame=frame,
        anchor_marker=anchor_marker,
        screenshot_path=args.screenshot,
    )
    viewer.run()


if __name__ == "__main__":
    main()
