# journal/aizu_demo/make_cells.py

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

# ------------------------------------------------------------
# Make repository root importable even without pip install -e .
# ------------------------------------------------------------
THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from aizu_sites import (
    SITES,
    AGL_LEVELS_M,
    XY_ZOOMS,
    VERTICAL_HALF_THICKNESS_M,
)

from bst_id.encoder import BSTIDEncoder
from bst_id.encoder_tile_calculator import (
    encode_x_tile_index,
    encode_y_tile_index,
    encode_f_tile_index,
)
from bst_id.constants import H_RANGE


# ============================================================
# Configuration
# ============================================================

# Current BST-ID altitude range is approximately 32767 m.
#
# z=11 -> about 16 m
# z=12 -> about  8 m
# z=13 -> about  4 m
#
# z=12 is a reasonable first choice for "about 10 m".
H_ZOOM = 12

OUTPUT_FILE = THIS_DIR / "aizu_cells.json"

EARTH_RADIUS_M = 6378137.0


# ============================================================
# Geographic conversion helpers
# ============================================================

def x_tile_bounds(tile_x: int, zoom: int) -> tuple[float, float]:
    """
    Return west/east longitude bounds of a Web-Mercator X tile.
    """
    n = 1 << zoom
    west = tile_x / n * 360.0 - 180.0
    east = (tile_x + 1) / n * 360.0 - 180.0
    return west, east


def tile_y_to_lat(tile_y: int, zoom: int) -> float:
    """
    Convert Web-Mercator tile Y boundary to latitude.
    """
    n = 1 << zoom
    value = math.pi * (1.0 - 2.0 * tile_y / n)
    return math.degrees(math.atan(math.sinh(value)))


def y_tile_bounds(tile_y: int, zoom: int) -> tuple[float, float]:
    """
    Return south/north latitude bounds of a Web-Mercator Y tile.

    Tile Y grows southward.
    """
    north = tile_y_to_lat(tile_y, zoom)
    south = tile_y_to_lat(tile_y + 1, zoom)
    return south, north


def h_tile_bounds(tile_h: int, zoom: int) -> tuple[float, float]:
    """
    Approximate altitude interval represented by an H-prefix.

    The current BST-ID implementation first maps height to a 32-bit
    altitude value and then retains the upper 'zoom' bits.
    """
    full_bits = 32
    shift = full_bits - zoom

    scale = (1 << full_bits) - 1
    h_min_domain, h_max_domain = H_RANGE
    span = h_max_domain - h_min_domain

    full_low = tile_h << shift
    full_high_exclusive = (tile_h + 1) << shift

    h_low = h_min_domain + (full_low / scale) * span
    h_high = h_min_domain + (full_high_exclusive / scale) * span

    return h_low, h_high


def haversine_m(
    lat1: float,
    lon1: float,
    lat2: float,
    lon2: float,
) -> float:
    """
    Approximate ground distance between two WGS84 lat/lon points.
    """
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    a = (
        math.sin(dphi / 2.0) ** 2
        + math.cos(phi1)
        * math.cos(phi2)
        * math.sin(dlambda / 2.0) ** 2
    )

    return 2.0 * EARTH_RADIUS_M * math.asin(math.sqrt(a))


# ============================================================
# BST-ID generation
# ============================================================

def make_cell(
    site_key: str,
    site: dict,
    agl_m: float,
    xy_zoom: int,
) -> dict:
    lat = float(site["lat"])
    lon = float(site["lon"])
    ground_m = float(site["ground_m"])

    # Center altitude of the intended airspace band.
    altitude_center_m = ground_m + agl_m

    # Desired Cesium visualization band.
    display_bottom_m = (
        altitude_center_m - VERTICAL_HALF_THICKNESS_M
    )
    display_top_m = (
        altitude_center_m + VERTICAL_HALF_THICKNESS_M
    )

    # Encode one representative XYH BST-ID at the center altitude.
    bst_id, bit_length = BSTIDEncoder.encode(
        x=lon,
        y=lat,
        f=altitude_center_m,
        t_unix=None,
        zoom_x=xy_zoom,
        zoom_y=xy_zoom,
        zoom_f=H_ZOOM,
        zoom_t=0,
    )

    # Recover tile indices explicitly for visualization.
    tile_x = encode_x_tile_index(
        lon,
        -180.0,
        180.0,
        xy_zoom,
    )

    tile_y = encode_y_tile_index(
        lat,
        xy_zoom,
    )

    tile_h = encode_f_tile_index(
        altitude_center_m,
        H_RANGE[0],
        H_RANGE[1],
        H_ZOOM,
    )

    west, east = x_tile_bounds(tile_x, xy_zoom)
    south, north = y_tile_bounds(tile_y, xy_zoom)
    bst_h_low, bst_h_high = h_tile_bounds(tile_h, H_ZOOM)

    # Approximate physical size of actual BST XY tile.
    center_lat = (south + north) / 2.0
    center_lon = (west + east) / 2.0

    width_m = haversine_m(
        center_lat,
        west,
        center_lat,
        east,
    )

    height_m = haversine_m(
        south,
        center_lon,
        north,
        center_lon,
    )

    return {
        "site_id": site_key,
        "site_name": site["name"],

        "source_point": {
            "latitude": lat,
            "longitude": lon,
            "ground_elevation_m": ground_m,
        },

        "operation_altitude": {
            "agl_center_m": agl_m,
            "amsl_center_m": altitude_center_m,

            # Requested visual slab: center +/- 10 m
            "display_bottom_m": display_bottom_m,
            "display_top_m": display_top_m,
            "display_thickness_m": (
                2.0 * VERTICAL_HALF_THICKNESS_M
            ),
        },

        "bst_id": {
            "decimal": str(bst_id),
            "hex": hex(bst_id),
            "binary": bin(bst_id),
            "bit_length": bit_length,

            "zoom": {
                "x": xy_zoom,
                "y": xy_zoom,
                "h": H_ZOOM,
            },

            "tile_index": {
                "x": tile_x,
                "y": tile_y,
                "h": tile_h,
            },
        },

        "bst_cell_bounds": {
            "west": west,
            "east": east,
            "south": south,
            "north": north,

            # Actual H-prefix interval of the representative ID.
            "h_bottom_m": bst_h_low,
            "h_top_m": bst_h_high,

            "approx_width_m": width_m,
            "approx_height_m": height_m,
            "approx_h_resolution_m": bst_h_high - bst_h_low,
        },

        # Ready to feed into Cesium RectangleGraphics.
        "cesium": {
            "rectangle": {
                "west": west,
                "south": south,
                "east": east,
                "north": north,
            },

            # Initially use the visually clearer requested +/-10 m band.
            "height": display_bottom_m,
            "extrudedHeight": display_top_m,
        },
    }


# ============================================================
# Main
# ============================================================

def main() -> None:
    cells = []

    ground_values = [
        float(site["ground_m"])
        for site in SITES.values()
    ]

    metadata = {
        "description": (
            "Aizuwakamatsu BST-ID / Cesium demonstration cells"
        ),

        "xy_zooms": XY_ZOOMS,
        "h_zoom": H_ZOOM,
        "agl_levels_m": AGL_LEVELS_M,

        "display_vertical_half_thickness_m":
            VERTICAL_HALF_THICKNESS_M,

        "ground_elevation": {
            "minimum_m": min(ground_values),
            "maximum_m": max(ground_values),
            "difference_m": max(ground_values) - min(ground_values),
        },

        "notes": [
            (
                "XY bounds are the actual grid-aligned BST-ID/Web-Mercator "
                "cell bounds, not a symmetric box around each source point."
            ),
            (
                "Cesium display height uses the requested AGL center "
                "+/- vertical half-thickness."
            ),
            (
                "The BST-ID H prefix is independently encoded at H zoom 12; "
                "its actual altitude interval is also included."
            ),
            (
                "Ground elevations are input values provided from "
                "GSI map elevation readings."
            ),
        ],
    }

    for site_key, site in SITES.items():
        for agl_m in AGL_LEVELS_M:
            for xy_zoom in XY_ZOOMS:
                cells.append(
                    make_cell(
                        site_key=site_key,
                        site=site,
                        agl_m=agl_m,
                        xy_zoom=xy_zoom,
                    )
                )

    output = {
        "metadata": metadata,
        "sites": SITES,
        "cells": cells,
    }

    OUTPUT_FILE.write_text(
        json.dumps(
            output,
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    print()
    print("BST-ID Aizu demonstration dataset")
    print("---------------------------------")
    print(f"Sites       : {len(SITES)}")
    print(f"AGL levels  : {AGL_LEVELS_M}")
    print(f"XY zooms    : {XY_ZOOMS}")
    print(f"H zoom      : {H_ZOOM}")
    print(f"Cells       : {len(cells)}")
    print()

    print(
        "Ground elevation range: "
        f"{metadata['ground_elevation']['minimum_m']:.1f} - "
        f"{metadata['ground_elevation']['maximum_m']:.1f} m "
        f"(difference "
        f"{metadata['ground_elevation']['difference_m']:.1f} m)"
    )

    print()
    print("Sites:")

    base_ground = min(ground_values)

    for site_key, site in SITES.items():
        delta = float(site["ground_m"]) - base_ground

        print(
            f"  {site['name']:28s} "
            f"{site['ground_m']:6.1f} m "
            f"(+{delta:4.1f} m from minimum)"
        )

    print()
    print(f"JSON written to:")
    print(f"  {OUTPUT_FILE}")
    print()

    # Show the z=20 / AGL=50 sample, useful for visual sanity check.
    print("Example: XY zoom 20, AGL 50 m")
    print("---------------------------------")

    for cell in cells:
        if (
            cell["bst_id"]["zoom"]["x"] == 20
            and cell["operation_altitude"]["agl_center_m"] == 50
        ):
            b = cell["bst_cell_bounds"]
            h = cell["operation_altitude"]

            print(
                f"{cell['site_name']}: "
                f"{b['approx_width_m']:.1f} x "
                f"{b['approx_height_m']:.1f} m, "
                f"AMSL "
                f"{h['display_bottom_m']:.1f}-"
                f"{h['display_top_m']:.1f} m"
            )


if __name__ == "__main__":
    main()