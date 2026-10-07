# journal/aizu_demo/generate_utm_scenario.py
from __future__ import annotations

import bisect
import json
import math
import random
import sys
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from aizu_sites import AGL_LEVELS_M, SITES
from bst_id.constants import H_RANGE, X_RANGE
from bst_id.encoder_tile_calculator import (
    encode_f_tile_index,
    encode_x_tile_index,
    encode_y_tile_index,
)
from bst_id.region_algebra import Cell

# =============================================================================
# Scenario configuration
# =============================================================================

JST = timezone(timedelta(hours=9))
SCENARIO_DATE = (2026, 10, 6)
SIM_START = datetime(*SCENARIO_DATE, 9, 0, 0, tzinfo=JST)
SIM_END = datetime(*SCENARIO_DATE, 11, 0, 0, tzinfo=JST)
LAST_DEPARTURE = datetime(*SCENARIO_DATE, 10, 40, 0, tzinfo=JST)

DEPARTURE_INTERVAL_S = 180
SAMPLE_INTERVAL_S = 1

CLIMB_SPEED_MPS = 3.0
DESCENT_SPEED_MPS = 3.0

# Three nominal cruise-speed classes.  Flights are assigned approximately
# evenly across the scenario, then shuffled deterministically with RANDOM_SEED.
CRUISE_SPEED_OPTIONS_MPS = (5.0, 8.0, 10.0)

# Speed-dependent route-following behavior.  Faster flights use a longer
# look-ahead distance and slightly larger residual disturbance; slower flights
# are calmer and track with a shorter spatial look-ahead.  This creates visibly
# different motion without allowing random-walk departure from the route.
CRUISE_PROFILES = {
    5.0: {
        "lookahead_m": 70.0,
        "heading_sigma_deg": 0.18,
        "heading_limit_deg": 1.1,
        "speed_sigma": 0.010,
        "speed_factor_min": 0.94,
        "speed_factor_max": 1.06,
        "lateral_gust_sigma_m": 0.020,
    },
    8.0: {
        "lookahead_m": 112.0,
        "heading_sigma_deg": 0.24,
        "heading_limit_deg": 1.4,
        "speed_sigma": 0.013,
        "speed_factor_min": 0.92,
        "speed_factor_max": 1.08,
        "lateral_gust_sigma_m": 0.025,
    },
    10.0: {
        "lookahead_m": 140.0,
        "heading_sigma_deg": 0.30,
        "heading_limit_deg": 1.8,
        "speed_sigma": 0.015,
        "speed_factor_min": 0.90,
        "speed_factor_max": 1.10,
        "lateral_gust_sigma_m": 0.030,
    },
}

BASE_XY_ZOOM = 20
MARGIN_XY_ZOOM = 21
H_ZOOM = 12

# Route-volume definition.  A 30 m route means exactly 20-40 m AGL,
# a 50 m route means 40-60 m AGL, and a 100 m route means 90-110 m AGL.
VERTICAL_BUBBLE_M = 10.0

# Each *nominal route sample* is valid from 5 min before to 5 min after
# the expected passage time.  This is deliberately NOT a whole-flight window.
TIME_BUBBLE_S = 300

# The route reservation covers the cruise corridor.  Takeoff/descent are shown
# by the moving UAS, but are not judged against the cruise corridor.  They can
# later be added as a separate departure/arrival operation volume if desired.
RESERVATION_PHASES = {"cruise"}

RANDOM_SEED = 20261006

# Actual-flight disturbance / guidance model.
# The aircraft is NOT allowed to random-walk away from the route.  A small
# correlated heading disturbance is superimposed on a look-ahead guidance
# command that continuously steers cross-track error back toward zero.
HEADING_NOISE_RHO = 0.82
TRACK_CORRECTION_LIMIT_DEG = 10.0
COMMAND_HEADING_LIMIT_DEG = 11.0

SPEED_RHO = 0.90
VERTICAL_SPEED_RHO = 0.85
VERTICAL_SPEED_INNOVATION_SIGMA = 0.05
VERTICAL_SPEED_FACTOR_MIN = 0.70
VERTICAL_SPEED_FACTOR_MAX = 1.30
CROSS_TRACK_LIMIT_M = 120.0
CROSS_TRACK_REVERSION = 0.998
CRUISE_VERTICAL_TRACK_RATE_MPS = 2.0
CRUISE_VERTICAL_NOISE_SIGMA_M = 0.15
TAKEOFF_LANDING_DRIFT_SIGMA_M = 0.08

EARTH_RADIUS_M = 6378137.0
OUTPUT_FILE = THIS_DIR / "utm_scenario.json"

Interval = Tuple[datetime, datetime]
Atom = Tuple[int, int, int]  # x21, y21, h12
XY = Tuple[int, int]

# =============================================================================
# General helpers
# =============================================================================

def clamp(value: float, low: float, high: float) -> float:
    return max(low, min(high, value))


def iso_jst(dt: datetime) -> str:
    return dt.astimezone(JST).isoformat(timespec="seconds")


def iso_utc(dt: datetime) -> str:
    return dt.astimezone(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def epoch_s(dt: datetime) -> int:
    return int(dt.timestamp())


def normalize_heading_deg(value: float) -> float:
    return value % 360.0


def haversine_m(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    p1 = math.radians(lat1)
    p2 = math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = (
        math.sin(dp / 2.0) ** 2
        + math.cos(p1) * math.cos(p2) * math.sin(dl / 2.0) ** 2
    )
    return 2.0 * EARTH_RADIUS_M * math.asin(math.sqrt(a))


def initial_bearing_deg(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    p1 = math.radians(lat1)
    p2 = math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return normalize_heading_deg(math.degrees(math.atan2(y, x)))


def destination_point(lat: float, lon: float, bearing_deg: float, distance_m: float) -> Tuple[float, float]:
    if distance_m == 0.0:
        return lat, lon
    brg = math.radians(bearing_deg)
    p1 = math.radians(lat)
    l1 = math.radians(lon)
    d = distance_m / EARTH_RADIUS_M
    p2 = math.asin(
        math.sin(p1) * math.cos(d)
        + math.cos(p1) * math.sin(d) * math.cos(brg)
    )
    l2 = l1 + math.atan2(
        math.sin(brg) * math.sin(d) * math.cos(p1),
        math.cos(d) - math.sin(p1) * math.sin(p2),
    )
    lon2 = (math.degrees(l2) + 540.0) % 360.0 - 180.0
    return math.degrees(p2), lon2


def interpolate_route_point(
    origin: Mapping[str, float],
    destination: Mapping[str, float],
    bearing_deg: float,
    route_distance_m: float,
    along_track_m: float,
) -> Tuple[float, float, float, float]:
    if route_distance_m <= 0.0:
        return float(origin["lat"]), float(origin["lon"]), float(origin["ground_m"]), 0.0
    fraction = clamp(along_track_m / route_distance_m, 0.0, 1.0)
    lat, lon = destination_point(
        float(origin["lat"]),
        float(origin["lon"]),
        bearing_deg,
        route_distance_m * fraction,
    )
    ground = (
        (1.0 - fraction) * float(origin["ground_m"])
        + fraction * float(destination["ground_m"])
    )
    return lat, lon, ground, fraction


def merge_intervals(intervals: Iterable[Interval]) -> List[Interval]:
    items = sorted((a, b) for a, b in intervals if b > a)
    if not items:
        return []
    out: List[List[datetime]] = [[items[0][0], items[0][1]]]
    for start, end in items[1:]:
        if start <= out[-1][1]:
            if end > out[-1][1]:
                out[-1][1] = end
        else:
            out.append([start, end])
    return [(a, b) for a, b in out]


def interval_contains(intervals: Sequence[Interval], t: datetime) -> bool:
    # The merged lists are normally tiny, so a short linear scan is cheaper and
    # clearer than building another index.
    for start, end in intervals:
        if start <= t < end:
            return True
        if start > t:
            break
    return False


def interval_json(intervals: Sequence[Interval]) -> List[Dict[str, object]]:
    return [
        {
            "start_jst": iso_jst(a),
            "end_jst": iso_jst(b),
            "start_utc": iso_utc(a),
            "end_utc": iso_utc(b),
            "start_epoch_s": epoch_s(a),
            "end_epoch_s": epoch_s(b),
        }
        for a, b in intervals
    ]

# =============================================================================
# BST-ID grid helpers
# =============================================================================

def x_tile_bounds(tile_x: int, zoom: int) -> Tuple[float, float]:
    n = 1 << zoom
    return tile_x / n * 360.0 - 180.0, (tile_x + 1) / n * 360.0 - 180.0


def tile_y_to_lat(tile_y: int, zoom: int) -> float:
    n = 1 << zoom
    value = math.pi * (1.0 - 2.0 * tile_y / n)
    return math.degrees(math.atan(math.sinh(value)))


def y_tile_bounds(tile_y: int, zoom: int) -> Tuple[float, float]:
    return tile_y_to_lat(tile_y + 1, zoom), tile_y_to_lat(tile_y, zoom)


def h_tile_bounds(tile_h: int, zoom: int) -> Tuple[float, float]:
    full_bits = 32
    shift = full_bits - zoom
    scale = (1 << full_bits) - 1
    h_min_domain, h_max_domain = H_RANGE
    span = h_max_domain - h_min_domain
    low = tile_h << shift
    high = (tile_h + 1) << shift
    return (
        h_min_domain + (low / scale) * span,
        h_min_domain + (high / scale) * span,
    )


def h_tiles_covering(low_m: float, high_m: float, zoom: int) -> range:
    lo = clamp(low_m, H_RANGE[0], H_RANGE[1])
    hi = clamp(high_m, H_RANGE[0], H_RANGE[1])
    if hi < lo:
        lo, hi = hi, lo
    first = encode_f_tile_index(lo, H_RANGE[0], H_RANGE[1], zoom)
    last = encode_f_tile_index(max(lo, hi - 1e-9), H_RANGE[0], H_RANGE[1], zoom)
    return range(first, last + 1)


def xy21_cells_for_base(ix20: int, iy20: int) -> Set[XY]:
    """z20 interior (4 z21 children) + 8 z21 N/E/S/W outward neighbors."""
    x0 = ix20 << 1
    y0 = iy20 << 1
    n = 1 << MARGIN_XY_ZOOM
    out: Set[XY] = set()
    for dx in (0, 1):
        for dy in (0, 1):
            out.add((x0 + dx, y0 + dy))
    for yy in (y0, y0 + 1):
        out.add(((x0 - 1) % n, yy))
        out.add(((x0 + 2) % n, yy))
    for xx in (x0, x0 + 1):
        if y0 - 1 >= 0:
            out.add((xx, y0 - 1))
        if y0 + 2 < n:
            out.add((xx, y0 + 2))
    return out


def reservation_query_atom(lon: float, lat: float, altitude_m: float) -> Atom:
    return (
        encode_x_tile_index(lon, X_RANGE[0], X_RANGE[1], MARGIN_XY_ZOOM),
        encode_y_tile_index(lat, MARGIN_XY_ZOOM),
        encode_f_tile_index(altitude_m, H_RANGE[0], H_RANGE[1], H_ZOOM),
    )


def atom_json(atom: Atom) -> Dict[str, object]:
    x21, y21, h12 = atom
    west, east = x_tile_bounds(x21, MARGIN_XY_ZOOM)
    south, north = y_tile_bounds(y21, MARGIN_XY_ZOOM)
    h_low, h_high = h_tile_bounds(h12, H_ZOOM)
    return {
        "x_zoom": MARGIN_XY_ZOOM,
        "y_zoom": MARGIN_XY_ZOOM,
        "h_zoom": H_ZOOM,
        "x_index": x21,
        "y_index": y21,
        "h_index": h12,
        "bounds": {
            "west": west,
            "east": east,
            "south": south,
            "north": north,
            "h_bottom_m": h_low,
            "h_top_m": h_high,
        },
    }

# =============================================================================
# Minimal dyadic prefix cover for Unix time
# =============================================================================

def _largest_power_of_two_leq(value: int) -> int:
    return 1 << (value.bit_length() - 1)


def time_prefix_cover(start_dt: datetime, end_dt: datetime) -> List[Dict[str, object]]:
    start = epoch_s(start_dt)
    end = epoch_s(end_dt)
    if end <= start:
        return []
    out = []
    cur = start
    while cur < end:
        remaining = end - cur
        aligned = (cur & -cur) if cur else (1 << 32)
        block = min(aligned, _largest_power_of_two_leq(remaining))
        while block > remaining:
            block >>= 1
        zoom = 32 - int(math.log2(block))
        if zoom < 1:
            zoom = 1
            block = 1 << 31
        index = cur >> (32 - zoom)
        out.append({
            "zoom_t": zoom,
            "t_index": index,
            "duration_s": block,
            "start_utc": iso_utc(datetime.fromtimestamp(cur, tz=timezone.utc)),
            "end_utc": iso_utc(datetime.fromtimestamp(cur + block, tz=timezone.utc)),
        })
        cur += block
    return out

# =============================================================================
# Nominal flight
# =============================================================================

def generate_nominal_trajectory(
    departure: datetime,
    origin: Mapping[str, float],
    destination: Mapping[str, float],
    target_agl_m: float,
    cruise_speed_mps: float,
) -> Tuple[List[Dict[str, object]], Dict[str, float]]:
    distance_m = haversine_m(
        float(origin["lat"]), float(origin["lon"]),
        float(destination["lat"]), float(destination["lon"]),
    )
    bearing_deg = initial_bearing_deg(
        float(origin["lat"]), float(origin["lon"]),
        float(destination["lat"]), float(destination["lon"]),
    )
    climb_duration = target_agl_m / CLIMB_SPEED_MPS
    cruise_duration = distance_m / cruise_speed_mps
    descent_duration = target_agl_m / DESCENT_SPEED_MPS
    total = climb_duration + cruise_duration + descent_duration
    nominal_duration_s = int(math.ceil(total))
    samples: List[Dict[str, object]] = []

    for sec in range(nominal_duration_s + 1):
        t = departure + timedelta(seconds=sec)
        if t > SIM_END:
            break
        e = float(sec)
        if e < climb_duration:
            phase = "climb"
            lat, lon = float(origin["lat"]), float(origin["lon"])
            ground = float(origin["ground_m"])
            agl = min(target_agl_m, e * CLIMB_SPEED_MPS)
            hs, vs = 0.0, CLIMB_SPEED_MPS
        elif e < climb_duration + cruise_duration:
            phase = "cruise"
            along = clamp((e - climb_duration) * cruise_speed_mps, 0.0, distance_m)
            lat, lon, ground, _ = interpolate_route_point(
                origin, destination, bearing_deg, distance_m, along
            )
            agl = target_agl_m
            hs = cruise_speed_mps
            vs = (
                float(destination["ground_m"]) - float(origin["ground_m"])
            ) / max(cruise_duration, 1e-9)
        elif e < total:
            phase = "descent"
            lat, lon = float(destination["lat"]), float(destination["lon"])
            ground = float(destination["ground_m"])
            agl = max(0.0, target_agl_m - (e - climb_duration - cruise_duration) * DESCENT_SPEED_MPS)
            hs, vs = 0.0, -DESCENT_SPEED_MPS
        else:
            phase = "landed"
            lat, lon = float(destination["lat"]), float(destination["lon"])
            ground = float(destination["ground_m"])
            agl, hs, vs = 0.0, 0.0, 0.0

        samples.append({
            "time_jst": iso_jst(t),
            "time_utc": iso_utc(t),
            "epoch_s": epoch_s(t),
            "elapsed_s": sec,
            "phase": phase,
            "lat": lat,
            "lon": lon,
            "ground_est_m": ground,
            "agl_m": agl,
            "amsl_m": ground + agl,
            "horizontal_speed_mps": hs,
            "vertical_speed_mps": vs,
            "heading_deg": bearing_deg,
        })

    return samples, {
        "distance_m": distance_m,
        "bearing_deg": bearing_deg,
        "climb_duration_s": climb_duration,
        "cruise_duration_s": cruise_duration,
        "descent_duration_s": descent_duration,
        "continuous_duration_s": total,
        "nominal_duration_s": nominal_duration_s,
        "cruise_speed_mps": cruise_speed_mps,
    }

# =============================================================================
# Dynamic route reservation
# =============================================================================

def build_dynamic_reservation(
    nominal_samples: Sequence[Dict[str, object]],
    target_agl_m: float,
) -> Tuple[
    Dict[Atom, List[Interval]],
    Dict[XY, List[Interval]],
    Dict[XY, Set[int]],
    List[Dict[str, object]],
    datetime,
    datetime,
]:
    """Build a moving 4-D reservation.

    Each cruise sample contributes one z20 route cell, expanded to z21 interior
    + N/E/S/W margin cells.  Only that sample's location is valid during
    [sample_time-5 min, sample_time+5 min).  Thus a current-time view guards
    the portion of the route expected within +/-5 min, not the whole route.

    The vertical route band is target AGL +/-10 m.  Thus:
       30 m route -> 20-40 m AGL
       50 m route -> 40-60 m AGL
      100 m route -> 90-110 m AGL
    BST-ID H cells are still encoded in absolute altitude using the locally
    estimated ground height.
    """
    raw_atom_intervals: Dict[Atom, List[Interval]] = defaultdict(list)
    raw_xy_intervals: Dict[XY, List[Interval]] = defaultdict(list)
    h_indices_by_xy: Dict[XY, Set[int]] = defaultdict(set)

    agl_low = max(0.0, target_agl_m - VERTICAL_BUBBLE_M)
    agl_high = target_agl_m + VERTICAL_BUBBLE_M

    for sample in nominal_samples:
        if str(sample["phase"]) not in RESERVATION_PHASES:
            continue

        t = datetime.fromtimestamp(int(sample["epoch_s"]), tz=timezone.utc).astimezone(JST)
        valid_from = max(SIM_START, t - timedelta(seconds=TIME_BUBBLE_S))
        valid_to = min(SIM_END, t + timedelta(seconds=TIME_BUBBLE_S))
        if valid_to <= valid_from:
            continue

        lon = float(sample["lon"])
        lat = float(sample["lat"])
        ground = float(sample["ground_est_m"])

        ix20 = encode_x_tile_index(lon, X_RANGE[0], X_RANGE[1], BASE_XY_ZOOM)
        iy20 = encode_y_tile_index(lat, BASE_XY_ZOOM)
        xy_cells = xy21_cells_for_base(ix20, iy20)

        absolute_low = ground + agl_low
        absolute_high = ground + agl_high
        h_indices = list(h_tiles_covering(absolute_low, absolute_high, H_ZOOM))

        for xy in xy_cells:
            raw_xy_intervals[xy].append((valid_from, valid_to))
            for ih12 in h_indices:
                atom = (xy[0], xy[1], ih12)
                raw_atom_intervals[atom].append((valid_from, valid_to))
                h_indices_by_xy[xy].add(ih12)

    atom_intervals = {k: merge_intervals(v) for k, v in raw_atom_intervals.items()}
    xy_intervals = {k: merge_intervals(v) for k, v in raw_xy_intervals.items()}

    all_intervals = [iv for vals in xy_intervals.values() for iv in vals]
    if all_intervals:
        reservation_start = min(a for a, _ in all_intervals)
        reservation_end = max(b for _, b in all_intervals)
    else:
        reservation_start = SIM_START
        reservation_end = SIM_START

    display_cells: List[Dict[str, object]] = []
    for (x21, y21), intervals in sorted(xy_intervals.items()):
        west, east = x_tile_bounds(x21, MARGIN_XY_ZOOM)
        south, north = y_tile_bounds(y21, MARGIN_XY_ZOOM)
        display_cells.append({
            "x21": x21,
            "y21": y21,
            "bounds": {
                "west": west,
                "east": east,
                "south": south,
                "north": north,
            },
            "agl_bottom_m": agl_low,
            "agl_top_m": agl_high,
            "h12_indices": sorted(h_indices_by_xy[(x21, y21)]),
            "intervals": interval_json(intervals),
        })

    return atom_intervals, xy_intervals, h_indices_by_xy, display_cells, reservation_start, reservation_end

# =============================================================================
# Actual trajectory + conformance
# =============================================================================

def generate_actual_trajectory(
    departure: datetime,
    origin: Mapping[str, float],
    destination: Mapping[str, float],
    target_agl_m: float,
    cruise_speed_mps: float,
    route_distance_m: float,
    route_bearing_deg: float,
    atom_intervals: Mapping[Atom, Sequence[Interval]],
    xy_intervals: Mapping[XY, Sequence[Interval]],
    reservation_start: datetime,
    reservation_end: datetime,
    rng: random.Random,
) -> Tuple[List[Dict[str, object]], Dict[str, int]]:
    if cruise_speed_mps not in CRUISE_PROFILES:
        raise ValueError(
            f"Unsupported cruise speed {cruise_speed_mps}; "
            f"expected one of {CRUISE_SPEED_OPTIONS_MPS}"
        )
    profile = CRUISE_PROFILES[cruise_speed_mps]
    track_lookahead_m = float(profile["lookahead_m"])
    heading_sigma_deg = float(profile["heading_sigma_deg"])
    heading_limit_deg = float(profile["heading_limit_deg"])
    speed_sigma = float(profile["speed_sigma"])
    speed_factor_min = float(profile["speed_factor_min"])
    speed_factor_max = float(profile["speed_factor_max"])
    lateral_gust_sigma_m = float(profile["lateral_gust_sigma_m"])

    state = "climb"
    current_time = departure
    altitude_m = float(origin["ground_m"])
    along_track_m = 0.0
    cross_track_m = 0.0
    heading_noise_deg = 0.0
    track_correction_deg = 0.0
    command_heading_offset_deg = 0.0
    speed_factor = 1.0
    vertical_speed_factor = 1.0

    samples: List[Dict[str, object]] = []
    stats = {
        "applicable": 0,
        "inside": 0,
        "outside": 0,
        "horizontal_violation": 0,
        "altitude_violation": 0,
        "time_violation": 0,
    }

    max_steps = int((SIM_END - departure).total_seconds()) + 1

    for step in range(max_steps):
        if state == "climb":
            center_lat, center_lon = float(origin["lat"]), float(origin["lon"])
            ground_est = float(origin["ground_m"])
            route_fraction = 0.0
            display_heading = route_bearing_deg
            horizontal_speed = 0.0
            vertical_speed = CLIMB_SPEED_MPS * vertical_speed_factor
        elif state == "cruise":
            center_lat, center_lon, ground_est, route_fraction = interpolate_route_point(
                origin, destination, route_bearing_deg, route_distance_m, along_track_m
            )
            display_heading = normalize_heading_deg(
                route_bearing_deg + command_heading_offset_deg
            )
            horizontal_speed = cruise_speed_mps * speed_factor
            vertical_speed = 0.0
        elif state == "descent":
            center_lat, center_lon = float(destination["lat"]), float(destination["lon"])
            ground_est = float(destination["ground_m"])
            route_fraction = 1.0
            display_heading = route_bearing_deg
            horizontal_speed = 0.0
            vertical_speed = -DESCENT_SPEED_MPS * vertical_speed_factor
        else:
            center_lat, center_lon = float(destination["lat"]), float(destination["lon"])
            ground_est = float(destination["ground_m"])
            route_fraction = 1.0
            display_heading = route_bearing_deg
            horizontal_speed = 0.0
            vertical_speed = 0.0

        lat, lon = destination_point(
            center_lat, center_lon, route_bearing_deg + 90.0, cross_track_m
        )

        applicable = state in RESERVATION_PHASES
        inside: Optional[bool]
        violation: Optional[str]
        query_atom = reservation_query_atom(lon, lat, altitude_m)

        if not applicable:
            inside = None
            violation = None
        else:
            stats["applicable"] += 1
            query_xy = (query_atom[0], query_atom[1])
            time_ok = reservation_start <= current_time < reservation_end
            xy_ok = interval_contains(xy_intervals.get(query_xy, ()), current_time)
            atom_ok = interval_contains(atom_intervals.get(query_atom, ()), current_time)

            if atom_ok:
                inside = True
                violation = None
                stats["inside"] += 1
            else:
                inside = False
                stats["outside"] += 1
                if not time_ok:
                    violation = "time"
                    stats["time_violation"] += 1
                elif not xy_ok:
                    violation = "horizontal"
                    stats["horizontal_violation"] += 1
                else:
                    violation = "altitude"
                    stats["altitude_violation"] += 1

        samples.append({
            "time_jst": iso_jst(current_time),
            "time_utc": iso_utc(current_time),
            "epoch_s": epoch_s(current_time),
            "elapsed_s": step,
            "phase": state,
            "lat": lat,
            "lon": lon,
            "ground_est_m": ground_est,
            "agl_m": altitude_m - ground_est,
            "amsl_m": altitude_m,
            "horizontal_speed_mps": horizontal_speed,
            "planned_cruise_speed_mps": cruise_speed_mps,
            "vertical_speed_mps": vertical_speed,
            "heading_deg": display_heading,
            # Kept for viewer/backward compatibility: this is the commanded
            # offset from route bearing (guidance correction + disturbance).
            "heading_error_deg": command_heading_offset_deg,
            "heading_disturbance_deg": heading_noise_deg,
            "track_correction_deg": track_correction_deg,
            "speed_factor": speed_factor,
            "cross_track_m": cross_track_m,
            "along_track_m": along_track_m,
            "route_fraction": route_fraction,
            "conformance_applicable": applicable,
            "inside_reservation": inside,
            "violation": violation,
            "reservation_atom": {
                "x21": query_atom[0],
                "y21": query_atom[1],
                "h12": query_atom[2],
            },
        })

        if current_time >= SIM_END or state == "landed":
            break

        speed_factor = 1.0 + SPEED_RHO * (speed_factor - 1.0) + rng.gauss(0.0, speed_sigma)
        speed_factor = clamp(speed_factor, speed_factor_min, speed_factor_max)
        vertical_speed_factor = (
            1.0
            + VERTICAL_SPEED_RHO * (vertical_speed_factor - 1.0)
            + rng.gauss(0.0, VERTICAL_SPEED_INNOVATION_SIGMA)
        )
        vertical_speed_factor = clamp(
            vertical_speed_factor, VERTICAL_SPEED_FACTOR_MIN, VERTICAL_SPEED_FACTOR_MAX
        )

        if state == "climb":
            cross_track_m += rng.gauss(0.0, TAKEOFF_LANDING_DRIFT_SIGMA_M)
            target_alt = float(origin["ground_m"]) + target_agl_m
            altitude_m = min(target_alt, altitude_m + CLIMB_SPEED_MPS * vertical_speed_factor)
            if altitude_m >= target_alt - 1e-6:
                altitude_m = target_alt
                heading_noise_deg = 0.0
                track_correction_deg = 0.0
                command_heading_offset_deg = 0.0
                state = "cruise"
        elif state == "cruise":
            # Small correlated disturbance (wind / controller residual).
            heading_noise_deg = (
                HEADING_NOISE_RHO * heading_noise_deg
                + rng.gauss(0.0, heading_sigma_deg)
            )
            heading_noise_deg = clamp(
                heading_noise_deg, -heading_limit_deg, heading_limit_deg
            )

            # Route-following guidance.  Positive cross-track is to the right
            # of the nominal course, so the correction is to the left.
            track_correction_deg = -math.degrees(
                math.atan2(cross_track_m, track_lookahead_m)
            )
            track_correction_deg = clamp(
                track_correction_deg,
                -TRACK_CORRECTION_LIMIT_DEG,
                TRACK_CORRECTION_LIMIT_DEG,
            )
            command_heading_offset_deg = clamp(
                track_correction_deg + heading_noise_deg,
                -COMMAND_HEADING_LIMIT_DEG,
                COMMAND_HEADING_LIMIT_DEG,
            )

            speed = cruise_speed_mps * speed_factor
            command_rad = math.radians(command_heading_offset_deg)
            along_track_m += speed * math.cos(command_rad)
            cross_track_m += speed * math.sin(command_rad)
            # Tiny lateral gust keeps the track non-ideal without producing a
            # visually alarming random walk.
            cross_track_m += rng.gauss(0.0, lateral_gust_sigma_m)
            cross_track_m = clamp(cross_track_m, -CROSS_TRACK_LIMIT_M, CROSS_TRACK_LIMIT_M)
            along_track_m = min(route_distance_m, along_track_m)

            _, _, target_ground, _ = interpolate_route_point(
                origin, destination, route_bearing_deg, route_distance_m, along_track_m
            )
            target_alt = target_ground + target_agl_m
            correction = clamp(
                target_alt - altitude_m,
                -CRUISE_VERTICAL_TRACK_RATE_MPS,
                CRUISE_VERTICAL_TRACK_RATE_MPS,
            )
            altitude_m += correction + rng.gauss(0.0, CRUISE_VERTICAL_NOISE_SIGMA_M)
            if along_track_m >= route_distance_m - 1e-6:
                along_track_m = route_distance_m
                state = "descent"
        elif state == "descent":
            cross_track_m *= 0.90
            cross_track_m += rng.gauss(0.0, TAKEOFF_LANDING_DRIFT_SIGMA_M)
            ground = float(destination["ground_m"])
            altitude_m = max(ground, altitude_m - DESCENT_SPEED_MPS * vertical_speed_factor)
            if altitude_m <= ground + 1e-6:
                altitude_m = ground
                cross_track_m = 0.0
                state = "landed"

        current_time += timedelta(seconds=SAMPLE_INTERVAL_S)

    return samples, stats

# =============================================================================
# Overlap
# =============================================================================

def sweep_overlap_segments(
    entries: Sequence[Tuple[datetime, datetime, str]]
) -> List[Tuple[datetime, datetime, Tuple[str, ...]]]:
    grouped: Dict[datetime, Dict[str, Set[str]]] = defaultdict(
        lambda: {"add": set(), "remove": set()}
    )
    for start, end, flight_id in entries:
        if end <= start:
            continue
        grouped[start]["add"].add(flight_id)
        grouped[end]["remove"].add(flight_id)

    active: Set[str] = set()
    previous: Optional[datetime] = None
    out: List[Tuple[datetime, datetime, Tuple[str, ...]]] = []
    for t in sorted(grouped):
        if previous is not None and t > previous and len(active) >= 2:
            out.append((previous, t, tuple(sorted(active))))
        for flight_id in grouped[t]["remove"]:
            active.discard(flight_id)
        for flight_id in grouped[t]["add"]:
            active.add(flight_id)
        previous = t
    return out


def max_concurrent_flight_windows(windows: Sequence[Tuple[datetime, datetime, str]]) -> int:
    events: List[Tuple[datetime, int]] = []
    for start, end, _ in windows:
        events.append((start, +1))
        events.append((end, -1))
    events.sort(key=lambda e: (e[0], e[1]))
    active = 0
    maximum = 0
    for _, delta in events:
        active += delta
        maximum = max(maximum, active)
    return maximum

# =============================================================================
# Main
# =============================================================================

def departure_slots() -> List[datetime]:
    out = []
    t = SIM_START
    while t <= LAST_DEPARTURE:
        out.append(t)
        t += timedelta(seconds=DEPARTURE_INTERVAL_S)
    return out


def balanced_cruise_speed_schedule(
    count: int,
    rng: random.Random,
) -> List[float]:
    """Return a deterministic, approximately balanced randomized speed list."""
    base = count // len(CRUISE_SPEED_OPTIONS_MPS)
    remainder = count % len(CRUISE_SPEED_OPTIONS_MPS)
    speeds: List[float] = []
    for i, speed in enumerate(CRUISE_SPEED_OPTIONS_MPS):
        speeds.extend([float(speed)] * (base + (1 if i < remainder else 0)))
    rng.shuffle(speeds)
    return speeds


def main() -> None:
    master_rng = random.Random(RANDOM_SEED)
    site_ids = list(SITES.keys())
    departures = departure_slots()
    cruise_speed_schedule = balanced_cruise_speed_schedule(len(departures), master_rng)

    flights_json: List[Dict[str, object]] = []
    global_atom_entries: Dict[Atom, List[Tuple[datetime, datetime, str]]] = defaultdict(list)
    flight_windows: List[Tuple[datetime, datetime, str]] = []

    total_nominal = 0
    total_actual = 0
    total_display_cells = 0
    total_dynamic_atoms = 0
    totals = defaultdict(int)

    for idx, departure in enumerate(departures, 1):
        flight_id = f"F{idx:03d}"
        origin_id = master_rng.choice(site_ids)
        destination_id = master_rng.choice([s for s in site_ids if s != origin_id])
        target_agl_m = float(master_rng.choice(AGL_LEVELS_M))
        cruise_speed_mps = float(cruise_speed_schedule[idx - 1])
        origin = SITES[origin_id]
        destination = SITES[destination_id]

        nominal, route_info = generate_nominal_trajectory(
            departure, origin, destination, target_agl_m, cruise_speed_mps
        )

        (
            atom_intervals,
            xy_intervals,
            h_indices_by_xy,
            display_cells,
            reservation_start,
            reservation_end,
        ) = build_dynamic_reservation(nominal, target_agl_m)

        if reservation_end > reservation_start:
            flight_windows.append((reservation_start, reservation_end, flight_id))

        for atom, intervals in atom_intervals.items():
            for start, end in intervals:
                global_atom_entries[atom].append((start, end, flight_id))

        actual_rng = random.Random(RANDOM_SEED + idx * 1009)
        actual, conformance = generate_actual_trajectory(
            departure=departure,
            origin=origin,
            destination=destination,
            target_agl_m=target_agl_m,
            cruise_speed_mps=cruise_speed_mps,
            route_distance_m=float(route_info["distance_m"]),
            route_bearing_deg=float(route_info["bearing_deg"]),
            atom_intervals=atom_intervals,
            xy_intervals=xy_intervals,
            reservation_start=reservation_start,
            reservation_end=reservation_end,
            rng=actual_rng,
        )

        total_nominal += len(nominal)
        total_actual += len(actual)
        total_display_cells += len(display_cells)
        total_dynamic_atoms += len(atom_intervals)
        for k, v in conformance.items():
            totals[k] += v

        # Compact inspection examples: one atom and one of its merged intervals.
        xyht_examples = []
        for atom, intervals in list(sorted(atom_intervals.items()))[:3]:
            for start, end in intervals[:1]:
                for tp in time_prefix_cover(start, end)[:3]:
                    x21, y21, h12 = atom
                    c = Cell.from_parts(
                        x=(MARGIN_XY_ZOOM, x21),
                        y=(MARGIN_XY_ZOOM, y21),
                        f=(H_ZOOM, h12),
                        t=(int(tp["zoom_t"]), int(tp["t_index"])),
                    )
                    xyht_examples.append({
                        "bst_id_xyht_hex": hex(c.to_id()),
                        "x21": x21,
                        "y21": y21,
                        "h12": h12,
                        "time_prefix": tp,
                    })

        flights_json.append({
            "flight_id": flight_id,
            "origin_id": origin_id,
            "origin_name": origin["name"],
            "destination_id": destination_id,
            "destination_name": destination["name"],
            "target_agl_m": target_agl_m,
            "planned_cruise_speed_mps": cruise_speed_mps,
            "speed_class": (
                "slow" if cruise_speed_mps == 5.0
                else "medium" if cruise_speed_mps == 8.0
                else "fast"
            ),
            "scheduled_departure_jst": iso_jst(departure),
            "scheduled_departure_utc": iso_utc(departure),
            "route": route_info,
            "reservation": {
                "start_jst": iso_jst(reservation_start),
                "end_jst": iso_jst(reservation_end),
                "start_utc": iso_utc(reservation_start),
                "end_utc": iso_utc(reservation_end),
                "time_bubble_before_s": TIME_BUBBLE_S,
                "time_bubble_after_s": TIME_BUBBLE_S,
                "vertical_bubble_m": VERTICAL_BUBBLE_M,
                "base_xy_zoom": BASE_XY_ZOOM,
                "margin_xy_zoom": MARGIN_XY_ZOOM,
                "h_zoom": H_ZOOM,
                "reserved_phases": sorted(RESERVATION_PHASES),
                "agl_bottom_m": max(0.0, target_agl_m - VERTICAL_BUBBLE_M),
                "agl_top_m": target_agl_m + VERTICAL_BUBBLE_M,
                "time_model": (
                    "Each nominal cruise sample reserves its local XY/H volume "
                    "from 5 min before to 5 min after the expected passage time; "
                    "the whole route is not activated as one mission-wide volume."
                ),
                "display_cells": display_cells,
                "display_cell_count": len(display_cells),
                "dynamic_atom_count": len(atom_intervals),
                # Exact time-varying XYH atoms for the Cesium operation-composition
                # view.  Intervals are stored as compact [start_epoch, end_epoch]
                # pairs to keep the JSON reasonably small.  The viewer performs:
                #   TimeSlice(t) -> Union -> Normalize(XYH) -> Project(H) -> Normalize(XY)
                # at the current Cesium time.
                "dynamic_atoms": [
                    {
                        "x21": atom[0],
                        "y21": atom[1],
                        "h12": atom[2],
                        "intervals": [
                            [epoch_s(start), epoch_s(end)]
                            for start, end in intervals
                        ],
                    }
                    for atom, intervals in sorted(atom_intervals.items())
                ],
                "xyht_examples": xyht_examples,
            },
            "nominal_trajectory": nominal,
            "actual_trajectory": actual,
            "conformance": conformance,
        })

        print(
            f"{flight_id} {origin['name']} -> {destination['name']} "
            f"AGL={int(target_agl_m)}m speed={cruise_speed_mps:g}m/s "
            f"depart={departure:%H:%M} "
            f"display={len(display_cells)} atoms={len(atom_intervals)} "
            f"actual={len(actual)} outside={conformance['outside']}"
        )

    overlap_segments_json: List[Dict[str, object]] = []
    overlap_atoms = 0
    max_overlap_coverage = 1

    for atom, entries in global_atom_entries.items():
        segments = sweep_overlap_segments(entries)
        if not segments:
            continue
        overlap_atoms += 1
        cell = atom_json(atom)
        for start, end, active_flights in segments:
            coverage = len(active_flights)
            max_overlap_coverage = max(max_overlap_coverage, coverage)
            overlap_segments_json.append({
                "cell": cell,
                "start_jst": iso_jst(start),
                "end_jst": iso_jst(end),
                "start_utc": iso_utc(start),
                "end_utc": iso_utc(end),
                "start_epoch_s": epoch_s(start),
                "end_epoch_s": epoch_s(end),
                "coverage": coverage,
                "flight_ids": list(active_flights),
            })

    applicable = totals["applicable"]
    inside_ratio = totals["inside"] / applicable if applicable else 0.0

    output = {
        "metadata": {
            "scenario": "BST-ID dynamic UTM reservation and conformance monitoring demo",
            "random_seed": RANDOM_SEED,
            "simulation_start_jst": iso_jst(SIM_START),
            "simulation_end_jst": iso_jst(SIM_END),
            "simulation_start_utc": iso_utc(SIM_START),
            "simulation_end_utc": iso_utc(SIM_END),
            "last_allowed_departure_jst": iso_jst(LAST_DEPARTURE),
            "departure_interval_s": DEPARTURE_INTERVAL_S,
            "flight_count": len(departures),
            "sample_interval_s": SAMPLE_INTERVAL_S,
            "speeds_mps": {
                "climb": CLIMB_SPEED_MPS,
                "cruise_options": list(CRUISE_SPEED_OPTIONS_MPS),
                "descent": DESCENT_SPEED_MPS,
            },
            "cruise_speed_assignment": {
                "method": "approximately balanced randomized assignment",
                "counts": {
                    str(int(speed)): cruise_speed_schedule.count(float(speed))
                    for speed in CRUISE_SPEED_OPTIONS_MPS
                },
            },
            "reservation": {
                "base_xy_zoom": BASE_XY_ZOOM,
                "margin_xy_zoom": MARGIN_XY_ZOOM,
                "h_zoom": H_ZOOM,
                "vertical_bubble_m": VERTICAL_BUBBLE_M,
                "time_bubble_before_s": TIME_BUBBLE_S,
                "time_bubble_after_s": TIME_BUBBLE_S,
                "reserved_phases": sorted(RESERVATION_PHASES),
            },
            "viewer_operation_composition": (
                "TimeSlice(t) -> Union(active XYH reservation atoms) -> "
                "Normalize(XYH by exact active-flight set and AGL band) -> "
                "Project(H to the operational AGL band) -> Normalize(XY for display)"
            ),
            "actual_flight_guidance": {
                "model": (
                    "speed-dependent route-following look-ahead guidance "
                    "with correlated heading disturbance"
                ),
                "track_correction_limit_deg": TRACK_CORRECTION_LIMIT_DEG,
                "profiles": {
                    str(int(speed)): CRUISE_PROFILES[float(speed)]
                    for speed in CRUISE_SPEED_OPTIONS_MPS
                },
            },
            "overlap_color_hint": {
                "coverage_1": "cyan",
                "coverage_2": "orange",
                "coverage_3_or_more": "red",
            },
            "conformance_color_hint": {
                "inside": "white drone",
                "outside": "red drone",
                "not_applicable": "gray/white drone during climb/descent",
            },
        },
        "sites": SITES,
        "summary": {
            "flight_count": len(departures),
            "cruise_speed_counts": {
                str(int(speed)): cruise_speed_schedule.count(float(speed))
                for speed in CRUISE_SPEED_OPTIONS_MPS
            },
            "total_nominal_samples": total_nominal,
            "total_actual_samples": total_actual,
            "total_display_cells_across_flights": total_display_cells,
            "total_dynamic_atoms_across_flights": total_dynamic_atoms,
            "unique_global_z21_h12_atoms": len(global_atom_entries),
            "overlap_atom_count": overlap_atoms,
            "overlap_segment_count": len(overlap_segments_json),
            "max_overlap_coverage": max_overlap_coverage,
            "max_concurrent_reservation_windows": max_concurrent_flight_windows(flight_windows),
            "conformance": {
                "applicable_samples": applicable,
                "inside_samples": totals["inside"],
                "outside_samples": totals["outside"],
                "inside_ratio": inside_ratio,
                "horizontal_violation_samples": totals["horizontal_violation"],
                "altitude_violation_samples": totals["altitude_violation"],
                "time_violation_samples": totals["time_violation"],
            },
        },
        "flights": flights_json,
        "overlap_segments": overlap_segments_json,
    }

    OUTPUT_FILE.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    print()
    print("Scenario complete (dynamic route reservation)")
    print("---------------------------------------------")
    print(f"Flights                    : {len(departures)}")
    print(f"Nominal samples            : {total_nominal}")
    print(f"Actual samples             : {total_actual}")
    print(f"Display cells              : {total_display_cells}")
    print(f"Dynamic z21/H12 atoms      : {total_dynamic_atoms}")
    print(f"Unique global z21/H12 atoms: {len(global_atom_entries)}")
    print(f"Overlap atoms              : {overlap_atoms}")
    print(f"Overlap segments           : {len(overlap_segments_json)}")
    print(f"Max overlap coverage       : {max_overlap_coverage}")
    print(f"Conformance applicable     : {applicable}")
    print(f"Conformance inside         : {totals['inside']}/{applicable} ({inside_ratio:.2%})")
    print(f"Output                     : {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
