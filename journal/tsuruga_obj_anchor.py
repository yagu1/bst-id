#!/usr/bin/env python3
"""
tsuruga_obj_anchor.py

Workflow:
  OBJ mesh
    -> uniform metric scale
    -> fixed yaw alignment
    -> surface point-cloud sampling
    -> pick ONE anchor point
    -> translate that point to local ENU origin
    -> save aligned PLY + transformation metadata

Coordinate convention after alignment:
    +X = East
    +Y = North
    +Z = Up
    units = metres

Default values are set for the current Tsuruga Castle experiment:
    scale   = 3.864
    yaw     = -13.64 deg
    anchor  = 37.4874741977344 N, 139.9299222365892 E

The geographic anchor is metadata only at this stage.  The point cloud remains
in a local metric ENU-like Cartesian frame, with the picked anchor mapped to
(0, 0, 0) by default.

Example:
    python tsuruga_obj_anchor.py tsuruga_20200216_260k.obj

Or explicitly:
    python tsuruga_obj_anchor.py \
        tsuruga_20200216_260k.obj \
        --scale 3.864 \
        --yaw-deg -13.64 \
        --anchor-lat 37.4874741977344 \
        --anchor-lon 139.9299222365892 \
        --n-points 300000 \
        --voxel 0.10 \
        --out tsuruga_enu_anchor.ply
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import open3d as o3d


def rotz(theta_rad: float) -> np.ndarray:
    """Rotation matrix about +Z. Positive angle is counter-clockwise in XY."""
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


def pick_one_point(pcd: o3d.geometry.PointCloud) -> int:
    """Interactively select exactly one point from a point cloud."""
    print()
    print("Open3D point picking")
    print("--------------------")
    print("Pick exactly ONE point to use as the geographic anchor.")
    print()
    print("  Shift + Left Click  : pick point")
    print("  Shift + Right Click : undo")
    print("  Q                    : close window")
    print()

    vis = o3d.visualization.VisualizerWithEditing()
    vis.create_window(window_name="Pick ONE anchor point")
    vis.add_geometry(pcd)
    vis.run()
    vis.destroy_window()

    picked = vis.get_picked_points()
    if len(picked) != 1:
        raise RuntimeError(
            f"Expected exactly 1 picked point, but got {len(picked)}."
        )
    return int(picked[0])


def make_frame(size: float = 10.0) -> o3d.geometry.TriangleMesh:
    # Open3D coordinate-frame colors:
    # X=red, Y=green, Z=blue.
    return o3d.geometry.TriangleMesh.create_coordinate_frame(
        size=size, origin=[0.0, 0.0, 0.0]
    )


def main() -> None:
    ap = argparse.ArgumentParser()

    ap.add_argument("obj", type=Path, help="Input OBJ mesh")

    ap.add_argument(
        "--out",
        type=Path,
        default=Path("tsuruga_enu_anchor.ply"),
        help="Output aligned point cloud",
    )
    ap.add_argument(
        "--transform-json",
        type=Path,
        default=Path("tsuruga_enu_anchor_transform.json"),
        help="Output transformation/anchor metadata",
    )

    # Current Tsuruga Castle values.
    ap.add_argument(
        "--scale",
        type=float,
        default=3.864,
        help="Uniform scale from OBJ units to metres",
    )
    ap.add_argument(
        "--yaw-deg",
        type=float,
        default=-13.64,
        help="Yaw about +Z in degrees; CCW is positive",
    )

    # Geographic anchor metadata.
    ap.add_argument(
        "--anchor-lat",
        type=float,
        default=37.4874741977344,
        help="Latitude of the picked anchor point [deg]",
    )
    ap.add_argument(
        "--anchor-lon",
        type=float,
        default=139.9299222365892,
        help="Longitude of the picked anchor point [deg]",
    )

    # Local target coordinate of the picked anchor.
    ap.add_argument(
        "--anchor-e",
        type=float,
        default=0.0,
        help="Target local East coordinate of anchor [m]",
    )
    ap.add_argument(
        "--anchor-n",
        type=float,
        default=0.0,
        help="Target local North coordinate of anchor [m]",
    )
    ap.add_argument(
        "--anchor-u",
        type=float,
        default=0.0,
        help="Target local Up coordinate of anchor [m]",
    )

    ap.add_argument(
        "--sample",
        choices=["uniform", "poisson"],
        default="uniform",
        help="Surface point sampling method",
    )
    ap.add_argument(
        "--n-points",
        type=int,
        default=300_000,
        help="Number of surface points sampled from mesh",
    )
    ap.add_argument(
        "--voxel",
        type=float,
        default=0.10,
        help="Optional point-cloud voxel downsample size [m]; 0 disables",
    )
    ap.add_argument(
        "--frame-size",
        type=float,
        default=10.0,
        help="Coordinate frame display size [m]",
    )
    ap.add_argument(
        "--no-show",
        action="store_true",
        help="Do not display final aligned point cloud",
    )

    args = ap.parse_args()

    # ------------------------------------------------------------------
    # 1. Load OBJ
    # ------------------------------------------------------------------
    mesh = o3d.io.read_triangle_mesh(
        str(args.obj),
        enable_post_processing=True,
    )
    if mesh.is_empty():
        raise RuntimeError(f"Could not load mesh: {args.obj}")

    mesh.compute_vertex_normals()

    bb0 = mesh.get_axis_aligned_bounding_box()
    print("Original mesh")
    print("-------------")
    print("vertices :", len(mesh.vertices))
    print("triangles:", len(mesh.triangles))
    print("min      :", np.asarray(bb0.min_bound))
    print("max      :", np.asarray(bb0.max_bound))
    print("extent   :", np.asarray(bb0.get_extent()))

    # ------------------------------------------------------------------
    # 2. Scale to metres
    # ------------------------------------------------------------------
    mesh.scale(args.scale, center=(0.0, 0.0, 0.0))

    # ------------------------------------------------------------------
    # 3. Apply fixed yaw
    # ------------------------------------------------------------------
    theta = math.radians(args.yaw_deg)
    R = rotz(theta)

    # Rotate mesh about the original OBJ origin.
    mesh.rotate(R, center=(0.0, 0.0, 0.0))

    bb1 = mesh.get_axis_aligned_bounding_box()
    print()
    print("After scale + yaw")
    print("-----------------")
    print(f"scale     : {args.scale:.9f}")
    print(f"yaw [deg] : {args.yaw_deg:.9f}")
    print("extent [m]:", np.asarray(bb1.get_extent()))

    # ------------------------------------------------------------------
    # 4. Sample mesh surface into point cloud
    # ------------------------------------------------------------------
    if args.sample == "poisson":
        pcd = mesh.sample_points_poisson_disk(
            number_of_points=args.n_points
        )
    else:
        pcd = mesh.sample_points_uniformly(
            number_of_points=args.n_points
        )

    if args.voxel > 0.0:
        pcd = pcd.voxel_down_sample(args.voxel)

    if pcd.is_empty():
        raise RuntimeError("Point cloud is empty after sampling.")

    # ------------------------------------------------------------------
    # 5. Pick ONE anchor point
    # ------------------------------------------------------------------
    anchor_idx = pick_one_point(pcd)
    pts = np.asarray(pcd.points)

    anchor_before_translation = pts[anchor_idx].copy()
    target_anchor = np.array(
        [args.anchor_e, args.anchor_n, args.anchor_u],
        dtype=float,
    )

    # Translation only: selected point -> requested local ENU coordinate.
    translation = target_anchor - anchor_before_translation

    pts_aligned = pts + translation
    pcd.points = o3d.utility.Vector3dVector(pts_aligned)

    # Selected anchor should now equal target_anchor (up to float precision).
    anchor_after_translation = np.asarray(pcd.points)[anchor_idx].copy()

    # ------------------------------------------------------------------
    # 6. Compute total transform from ORIGINAL OBJ coordinate to local ENU
    #
    # p_local = scale * R @ p_obj + translation
    # ------------------------------------------------------------------
    A = np.eye(4, dtype=float)
    A[:3, :3] = args.scale * R
    A[:3, 3] = translation

    # ------------------------------------------------------------------
    # 7. Report final bounds
    # ------------------------------------------------------------------
    bb2 = pcd.get_axis_aligned_bounding_box()

    print()
    print("Anchor")
    print("------")
    print("point index                  :", anchor_idx)
    print("anchor before translation [m]:", anchor_before_translation)
    print("translation [m]              :", translation)
    print("anchor after translation [m] :", anchor_after_translation)
    print("geographic anchor            :")
    print("  latitude  [deg]:", args.anchor_lat)
    print("  longitude [deg]:", args.anchor_lon)

    print()
    print("Final local point cloud")
    print("-----------------------")
    print("points     :", len(pcd.points))
    print("min [m]    :", np.asarray(bb2.min_bound))
    print("max [m]    :", np.asarray(bb2.max_bound))
    print("extent [m] :", np.asarray(bb2.get_extent()))
    print("axes       : +X East, +Y North, +Z Up")

    # ------------------------------------------------------------------
    # 8. Save aligned PLY
    # ------------------------------------------------------------------
    if not o3d.io.write_point_cloud(
        str(args.out),
        pcd,
        write_ascii=False,
    ):
        raise RuntimeError(f"Could not write point cloud: {args.out}")

    # ------------------------------------------------------------------
    # 9. Save metadata
    # ------------------------------------------------------------------
    metadata = {
        "input_obj": str(args.obj),
        "output_point_cloud": str(args.out),

        "source_to_metric": {
            "uniform_scale_m_per_obj_unit": args.scale,
            "yaw_deg_about_up": args.yaw_deg,
        },

        "coordinate_convention": {
            "x": "East",
            "y": "North",
            "z": "Up",
            "units": "metres",
        },

        "geographic_anchor": {
            "latitude_deg": args.anchor_lat,
            "longitude_deg": args.anchor_lon,
            "absolute_altitude_m": None,
            "note": (
                "Absolute altitude is intentionally unspecified at this stage. "
                "The picked physical anchor is mapped to local U=anchor_u."
            ),
        },

        "local_anchor_target_m": {
            "east": args.anchor_e,
            "north": args.anchor_n,
            "up": args.anchor_u,
        },

        "picked_anchor": {
            "point_index": anchor_idx,
            "position_after_scale_and_yaw_before_translation_m":
                anchor_before_translation.tolist(),
            "position_after_translation_m":
                anchor_after_translation.tolist(),
        },

        "translation_after_scale_and_yaw_m": translation.tolist(),

        "affine_transform_original_OBJ_to_local_ENU": A.tolist(),

        "sampling": {
            "method": args.sample,
            "requested_points": args.n_points,
            "voxel_downsample_m": args.voxel,
            "output_points": len(pcd.points),
        },

        "final_bounds_m": {
            "min": np.asarray(bb2.min_bound).tolist(),
            "max": np.asarray(bb2.max_bound).tolist(),
            "extent": np.asarray(bb2.get_extent()).tolist(),
        },
    }

    args.transform_json.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print()
    print("Saved")
    print("-----")
    print(args.out)
    print(args.transform_json)

    # ------------------------------------------------------------------
    # 10. Show result
    # ------------------------------------------------------------------
    if not args.no_show:
        frame = make_frame(args.frame_size)
        o3d.visualization.draw_geometries(
            [pcd, frame],
            window_name=(
                "Tsuruga Castle local frame "
                "(+X East, +Y North, +Z Up)"
            ),
        )


if __name__ == "__main__":
    main()
