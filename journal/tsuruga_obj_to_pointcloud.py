#!/usr/bin/env python3
"""
tsuruga_obj_to_pointcloud.py

Load an OBJ mesh with Open3D, sample it into a point cloud, apply a metric
scale, and optionally align the model to a local ENU-like frame:

    +X = East
    +Y = North
    +Z = Up

Two alignment modes are available:

1) Manual yaw:
       --yaw-deg <degrees>

2) Interactive north alignment:
       --pick-north

   Pick TWO points in the Open3D window:
     1. anchor point
     2. a point that should lie NORTH of the anchor

   The script translates the anchor to the requested local anchor coordinate
   (default: 0,0,0) and rotates about +Z so that the anchor->north vector
   points along +Y.

Example:
    python tsuruga_obj_to_pointcloud.py \
        tsuruga_20200216_260k.obj \
        --scale 3.864 \
        --n-points 300000 \
        --voxel 0.10 \
        --pick-north \
        --out tsuruga_aligned.ply
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import open3d as o3d


def rotz(theta_rad: float) -> np.ndarray:
    c = math.cos(theta_rad)
    s = math.sin(theta_rad)
    return np.array(
        [
            [c, -s, 0.0],
            [s,  c, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )


def pick_two_points(pcd: o3d.geometry.PointCloud) -> tuple[int, int]:
    print()
    print("Open3D point picking:")
    print("  Shift + Left Click : pick point")
    print("  Shift + Right Click: undo")
    print("  Q                  : close window")
    print()
    print("Pick exactly TWO points:")
    print("  1) anchor")
    print("  2) point that should be NORTH of anchor")
    print()

    vis = o3d.visualization.VisualizerWithEditing()
    vis.create_window(window_name="Pick anchor, then north-reference point")
    vis.add_geometry(pcd)
    vis.run()
    vis.destroy_window()

    picked = vis.get_picked_points()
    if len(picked) != 2:
        raise RuntimeError(
            f"Expected exactly 2 picked points, but got {len(picked)}."
        )
    return int(picked[0]), int(picked[1])


def make_frame(size: float = 10.0) -> o3d.geometry.TriangleMesh:
    # Open3D frame convention: X=red, Y=green, Z=blue.
    return o3d.geometry.TriangleMesh.create_coordinate_frame(
        size=size, origin=[0.0, 0.0, 0.0]
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("obj", type=Path, help="Input OBJ file")
    ap.add_argument("--out", type=Path, default=Path("tsuruga_aligned.ply"))
    ap.add_argument(
        "--transform-json",
        type=Path,
        default=Path("tsuruga_transform.json"),
    )

    ap.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="Uniform metric scale multiplier applied to OBJ coordinates",
    )
    ap.add_argument(
        "--yaw-deg",
        type=float,
        default=0.0,
        help="Manual yaw rotation about +Z, degrees (CCW positive)",
    )
    ap.add_argument(
        "--pick-north",
        action="store_true",
        help="Interactively pick anchor + north-reference point",
    )

    ap.add_argument(
        "--anchor-e",
        type=float,
        default=0.0,
        help="Target East coordinate for picked anchor [m]",
    )
    ap.add_argument(
        "--anchor-n",
        type=float,
        default=0.0,
        help="Target North coordinate for picked anchor [m]",
    )
    ap.add_argument(
        "--anchor-u",
        type=float,
        default=0.0,
        help="Target Up coordinate for picked anchor [m]",
    )

    ap.add_argument(
        "--sample",
        choices=["uniform", "poisson"],
        default="uniform",
    )
    ap.add_argument("--n-points", type=int, default=300_000)
    ap.add_argument(
        "--voxel",
        type=float,
        default=0.0,
        help="Optional voxel downsample size in metres; 0 disables",
    )
    ap.add_argument(
        "--frame-size",
        type=float,
        default=10.0,
        help="Displayed coordinate-frame size [m]",
    )
    ap.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open the final visualization window",
    )
    args = ap.parse_args()

    # ------------------------------------------------------------------
    # 1. Load mesh
    # ------------------------------------------------------------------
    mesh = o3d.io.read_triangle_mesh(str(args.obj), enable_post_processing=True)
    if mesh.is_empty():
        raise RuntimeError(f"Could not load mesh: {args.obj}")

    mesh.compute_vertex_normals()

    verts0 = np.asarray(mesh.vertices).copy()
    bb0 = mesh.get_axis_aligned_bounding_box()
    print("Original mesh:")
    print("  vertices :", len(mesh.vertices))
    print("  triangles:", len(mesh.triangles))
    print("  min      :", np.asarray(bb0.min_bound))
    print("  max      :", np.asarray(bb0.max_bound))
    print("  extent   :", np.asarray(bb0.get_extent()))

    # ------------------------------------------------------------------
    # 2. Uniform scale
    # ------------------------------------------------------------------
    mesh.scale(args.scale, center=(0.0, 0.0, 0.0))

    bb1 = mesh.get_axis_aligned_bounding_box()
    print()
    print(f"After scale x {args.scale:.9f}:")
    print("  extent [m]:", np.asarray(bb1.get_extent()))

    # ------------------------------------------------------------------
    # 3. Sample surface to point cloud
    # ------------------------------------------------------------------
    if args.sample == "poisson":
        pcd = mesh.sample_points_poisson_disk(number_of_points=args.n_points)
    else:
        pcd = mesh.sample_points_uniformly(number_of_points=args.n_points)

    if args.voxel > 0.0:
        pcd = pcd.voxel_down_sample(args.voxel)

    pts = np.asarray(pcd.points)
    if len(pts) == 0:
        raise RuntimeError("Point cloud is empty after sampling.")

    # Transformation bookkeeping.
    # x_final = R_total @ x_scaled + t_total
    R_total = np.eye(3)
    t_total = np.zeros(3)

    # ------------------------------------------------------------------
    # 4. Optional manual yaw
    # ------------------------------------------------------------------
    if abs(args.yaw_deg) > 1e-12:
        theta = math.radians(args.yaw_deg)
        Rm = rotz(theta)

        pts = (Rm @ pts.T).T
        pcd.points = o3d.utility.Vector3dVector(pts)

        R_total = Rm @ R_total
        t_total = Rm @ t_total

        print()
        print(f"Applied manual yaw: {args.yaw_deg:.6f} deg")

    # ------------------------------------------------------------------
    # 5. Interactive anchor + north alignment
    # ------------------------------------------------------------------
    picked_info = None

    if args.pick_north:
        i_anchor, i_north = pick_two_points(pcd)
        pts = np.asarray(pcd.points)

        anchor = pts[i_anchor].copy()
        north_pt = pts[i_north].copy()

        v = north_pt - anchor
        v[2] = 0.0
        norm_xy = np.linalg.norm(v[:2])
        if norm_xy < 1e-9:
            raise RuntimeError(
                "Picked anchor and north-reference have almost identical XY."
            )

        # Current XY azimuth from +X, target direction is +Y (= pi/2).
        phi = math.atan2(v[1], v[0])
        theta = (math.pi / 2.0) - phi
        Rn = rotz(theta)

        target_anchor = np.array(
            [args.anchor_e, args.anchor_n, args.anchor_u], dtype=float
        )

        # Rotate about the picked anchor and place the anchor at target coordinate.
        pts_aligned = (Rn @ (pts - anchor).T).T + target_anchor
        pcd.points = o3d.utility.Vector3dVector(pts_aligned)

        # Compose affine transform:
        # y = Rn*x + (target - Rn*anchor)
        t_step = target_anchor - Rn @ anchor
        R_total = Rn @ R_total
        t_total = Rn @ t_total + t_step

        picked_info = {
            "anchor_index": i_anchor,
            "north_reference_index": i_north,
            "anchor_before_alignment": anchor.tolist(),
            "north_reference_before_alignment": north_pt.tolist(),
            "computed_yaw_deg": math.degrees(theta),
            "target_anchor_ENU_m": target_anchor.tolist(),
        }

        print()
        print("Interactive alignment:")
        print("  anchor index       :", i_anchor)
        print("  north-ref index    :", i_north)
        print("  computed yaw [deg] :", math.degrees(theta))
        print("  anchor target [m]  :", target_anchor)

    # ------------------------------------------------------------------
    # 6. Report aligned bounds
    # ------------------------------------------------------------------
    bb2 = pcd.get_axis_aligned_bounding_box()
    print()
    print("Final point cloud:")
    print("  points      :", len(pcd.points))
    print("  min [m]     :", np.asarray(bb2.min_bound))
    print("  max [m]     :", np.asarray(bb2.max_bound))
    print("  extent [m]  :", np.asarray(bb2.get_extent()))
    print("  convention  : +X East, +Y North, +Z Up")

    # ------------------------------------------------------------------
    # 7. Save PLY + transform metadata
    # ------------------------------------------------------------------
    ok = o3d.io.write_point_cloud(str(args.out), pcd, write_ascii=False)
    if not ok:
        raise RuntimeError(f"Could not write point cloud: {args.out}")

    transform4 = np.eye(4)
    transform4[:3, :3] = R_total * args.scale
    transform4[:3, 3] = t_total

    meta = {
        "input_obj": str(args.obj),
        "output_point_cloud": str(args.out),
        "uniform_scale": args.scale,
        "manual_yaw_deg": args.yaw_deg,
        "pick_north": args.pick_north,
        "sample_method": args.sample,
        "requested_sample_points": args.n_points,
        "voxel_downsample_m": args.voxel,
        "coordinate_convention": {
            "x": "East",
            "y": "North",
            "z": "Up",
            "units": "metres",
        },
        "picked_alignment": picked_info,
        "affine_transform_from_original_OBJ": transform4.tolist(),
        "final_bounds_m": {
            "min": np.asarray(bb2.min_bound).tolist(),
            "max": np.asarray(bb2.max_bound).tolist(),
            "extent": np.asarray(bb2.get_extent()).tolist(),
        },
    }

    args.transform_json.write_text(
        json.dumps(meta, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print()
    print("Saved:")
    print(" ", args.out)
    print(" ", args.transform_json)

    # ------------------------------------------------------------------
    # 8. Visualize
    # ------------------------------------------------------------------
    if not args.no_show:
        frame = make_frame(args.frame_size)
        o3d.visualization.draw_geometries(
            [pcd, frame],
            window_name="Tsuruga Castle point cloud (+X East, +Y North, +Z Up)",
        )


if __name__ == "__main__":
    main()
