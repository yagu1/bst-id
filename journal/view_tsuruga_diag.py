#!/usr/bin/env python3
"""
view_tsuruga_diag.py

Simple Open3D viewer for the diagnostic PLY files produced by
check_tsuruga_1m_pipeline.py.

Keys:
  1 : classified overlay (all diagnostic categories)
  2 : raw shell only
  3 : closing-added cells only
  4 : fill-added cells only
  5 : fill-created exposed boundary only
  A : show all individual layers together
  C : clear all diagnostic layers
  R : reset camera
  P : save screenshot to tsuruga_diag_view.png
  Q / Esc : close

Color legend:
  gray   = raw shell
  orange = added by closing
  blue   = added by directional fill
  red    = added by directional fill AND exposed on final exterior
"""

from __future__ import annotations

import argparse
from pathlib import Path
import open3d as o3d


def load_cloud(path: Path) -> o3d.geometry.PointCloud:
    pcd = o3d.io.read_point_cloud(str(path))
    if pcd.is_empty():
        raise RuntimeError(f"Could not read {path}")
    return pcd


class Viewer:
    def __init__(self, args):
        self.args = args
        self.vis = o3d.visualization.VisualizerWithKeyCallback()

        self.layers = {
            "classified": load_cloud(args.classified),
            "raw": load_cloud(args.raw),
            "closed": load_cloud(args.closed),
            "fill": load_cloud(args.fill),
            "exposed": load_cloud(args.exposed),
        }

        self.visible = set()
        self.initialized = False

    def _add(self, name, reset=False):
        if name in self.visible:
            return
        geom = self.layers[name]
        try:
            self.vis.add_geometry(geom, reset_bounding_box=reset)
        except TypeError:
            self.vis.add_geometry(geom)
        self.visible.add(name)

    def _remove(self, name):
        if name not in self.visible:
            return
        geom = self.layers[name]
        try:
            self.vis.remove_geometry(geom, reset_bounding_box=False)
        except TypeError:
            self.vis.remove_geometry(geom)
        self.visible.discard(name)

    def _clear(self):
        for name in list(self.visible):
            self._remove(name)

    def show_only(self, name):
        self._clear()
        self._add(name, reset=True)

    def show_all_individual(self):
        self._clear()
        self._add("raw", reset=True)
        self._add("closed")
        self._add("fill")
        self._add("exposed")

    def cb1(self, vis):
        self.show_only("classified")
        return False

    def cb2(self, vis):
        self.show_only("raw")
        return False

    def cb3(self, vis):
        self.show_only("closed")
        return False

    def cb4(self, vis):
        self.show_only("fill")
        return False

    def cb5(self, vis):
        self.show_only("exposed")
        return False

    def cba(self, vis):
        self.show_all_individual()
        return False

    def cbc(self, vis):
        self._clear()
        return False

    def cbr(self, vis):
        ctr = self.vis.get_view_control()
        ctr.set_zoom(0.55)
        ctr.set_front([0.65, -0.55, -0.52])
        ctr.set_lookat([0.0, 0.0, 8.0])
        ctr.set_up([0.0, 0.0, 1.0])
        return False

    def cbp(self, vis):
        self.vis.capture_screen_image(str(self.args.screenshot), do_render=True)
        print("Saved screenshot:", self.args.screenshot)
        return False

    def run(self):
        self.vis.create_window(
            window_name="Tsuruga BST-ID diagnostic viewer",
            width=1400,
            height=900,
        )

        # Default: classified overlay.
        self._add("classified", reset=True)

        self.vis.register_key_callback(ord("1"), self.cb1)
        self.vis.register_key_callback(ord("2"), self.cb2)
        self.vis.register_key_callback(ord("3"), self.cb3)
        self.vis.register_key_callback(ord("4"), self.cb4)
        self.vis.register_key_callback(ord("5"), self.cb5)
        self.vis.register_key_callback(ord("A"), self.cba)
        self.vis.register_key_callback(ord("C"), self.cbc)
        self.vis.register_key_callback(ord("R"), self.cbr)
        self.vis.register_key_callback(ord("P"), self.cbp)

        opt = self.vis.get_render_option()
        opt.point_size = float(self.args.point_size)
        opt.background_color = [1.0, 1.0, 1.0]

        print()
        print("Tsuruga diagnostic viewer")
        print("=========================")
        print("1 : classified overlay")
        print("2 : raw shell only")
        print("3 : closing-added only")
        print("4 : fill-added only")
        print("5 : fill-created exposed boundary only")
        print("A : all individual layers together")
        print("C : clear")
        print("R : reset camera")
        print("P : screenshot")
        print("Q / Esc : close")
        print()
        print("Legend:")
        print("  gray   raw shell")
        print("  orange closing-added")
        print("  blue   fill-added")
        print("  red    fill-created exterior boundary")

        self.vis.run()
        self.vis.destroy_window()


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--classified",
        type=Path,
        default=Path("tsuruga_diag_classified.ply"),
    )
    ap.add_argument(
        "--raw",
        type=Path,
        default=Path("tsuruga_diag_raw_shell.ply"),
    )
    ap.add_argument(
        "--closed",
        type=Path,
        default=Path("tsuruga_diag_closed_added.ply"),
    )
    ap.add_argument(
        "--fill",
        type=Path,
        default=Path("tsuruga_diag_fill_added.ply"),
    )
    ap.add_argument(
        "--exposed",
        type=Path,
        default=Path("tsuruga_diag_fill_exposed_boundary.ply"),
    )
    ap.add_argument(
        "--point-size",
        type=float,
        default=5.0,
    )
    ap.add_argument(
        "--screenshot",
        type=Path,
        default=Path("tsuruga_diag_view.png"),
    )

    args = ap.parse_args()
    Viewer(args).run()


if __name__ == "__main__":
    main()
