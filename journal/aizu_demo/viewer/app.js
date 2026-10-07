// journal/aizu_demo/viewer/app.js
// Viewer for generate_utm_scenario.py (dynamic per-location +/-5 min reservation)

if (!window.CESIUM_ION_TOKEN) {
  throw new Error("CESIUM_ION_TOKEN is not defined. Check config.local.js");
}
Cesium.Ion.defaultAccessToken = window.CESIUM_ION_TOKEN;

const viewer = new Cesium.Viewer("cesiumContainer", {
  terrain: Cesium.Terrain.fromWorldTerrain(),
  animation: true,
  timeline: true,
  sceneModePicker: true,
  baseLayerPicker: false,
  geocoder: false,
  homeButton: false,
});

viewer.scene.globe.depthTestAgainstTerrain = true;

// -----------------------------------------------------------------------------
// Explicit GSI base map
// -----------------------------------------------------------------------------
//
// Do not rely on Cesium's default imagery provider.  The previous integrated
// build could render the terrain/objects while the default imagery failed,
// which makes the globe appear almost black.  Use GSI tiles explicitly.

viewer.scene.globe.baseColor = Cesium.Color.WHITE;

function addGsiBaseMap(kind = "std") {
  viewer.imageryLayers.removeAll();

  const isAerial = kind === "ort";
  const url = isAerial
    ? "https://cyberjapandata.gsi.go.jp/xyz/ort/{z}/{x}/{y}.jpg"
    : "https://cyberjapandata.gsi.go.jp/xyz/std/{z}/{x}/{y}.png";

  const layer = viewer.imageryLayers.addImageryProvider(
    new Cesium.UrlTemplateImageryProvider({
      url,
      credit: new Cesium.Credit("地理院タイル"),
      maximumLevel: 18,
    })
  );

  layer.alpha = 1.0;
  layer.show = true;
  return layer;
}

let gsiBaseLayer = addGsiBaseMap("std");

viewer.clock.clockRange = Cesium.ClockRange.CLAMPED;
viewer.clock.clockStep = Cesium.ClockStep.SYSTEM_CLOCK_MULTIPLIER;
viewer.clock.multiplier = 30;
viewer.clock.shouldAnimate = false;

const RELATIVE_TO_TERRAIN =
  Cesium.HeightReference.RELATIVE_TO_TERRAIN ??
  Cesium.HeightReference.RELATIVE_TO_GROUND;
const CLAMP_TO_TERRAIN =
  Cesium.HeightReference.CLAMP_TO_TERRAIN ??
  Cesium.HeightReference.CLAMP_TO_GROUND;

// -----------------------------------------------------------------------------
// Toolbar
// -----------------------------------------------------------------------------
let toolbar = document.getElementById("toolbar");
if (!toolbar) {
  toolbar = document.createElement("div");
  toolbar.id = "toolbar";
  document.body.appendChild(toolbar);
}

toolbar.innerHTML = `
  <div class="title">BST-ID UTM Reservation Demo</div>

  <div class="toolbar-controls">
    <label>
      Speed
      <select id="speedSelect">
        <option value="1">1×</option>
        <option value="10">10×</option>
        <option value="30" selected>30×</option>
        <option value="60">60×</option>
        <option value="120">120×</option>
      </select>
    </label>

    <label>
      Base map
      <select id="baseMapSelect">
        <option value="std" selected>GSI standard</option>
        <option value="ort">GSI aerial</option>
      </select>
    </label>

    <label><input id="showReservations" type="checkbox" checked> Reservations</label>
    <label><input id="overlapOnly" type="checkbox"> Overlap only</label>
    <label><input id="showLabels" type="checkbox" checked> Flight IDs</label>
    <label><input id="showRoutes" type="checkbox" checked> Routes</label>
    <label><input id="showBuildings" type="checkbox" checked> OSM 3D</label>

    <button id="homeView" type="button">Aizu view</button>
  </div>

  <div id="scenarioStatus" class="status-box">
    Loading scenario...
  </div>
`;

const showReservationsEl = document.getElementById("showReservations");
const overlapOnlyEl = document.getElementById("overlapOnly");
const showLabelsEl = document.getElementById("showLabels");
const showRoutesEl = document.getElementById("showRoutes");
const showBuildingsEl = document.getElementById("showBuildings");
const speedSelectEl = document.getElementById("speedSelect");
const baseMapSelectEl = document.getElementById("baseMapSelect");
const scenarioStatusEl = document.getElementById("scenarioStatus");

// -----------------------------------------------------------------------------
// Load scenario
// -----------------------------------------------------------------------------
const response = await fetch("../utm_scenario.json");
if (!response.ok) {
  throw new Error(`Failed to load utm_scenario.json: HTTP ${response.status}`);
}
const dataset = await response.json();

const simStart = Cesium.JulianDate.fromIso8601(dataset.metadata.simulation_start_utc);
const simStop = Cesium.JulianDate.fromIso8601(dataset.metadata.simulation_end_utc);
viewer.clock.startTime = simStart.clone();
viewer.clock.stopTime = simStop.clone();
viewer.clock.currentTime = simStart.clone();
viewer.timeline.zoomTo(simStart, simStop);

const jstFormatter = new Intl.DateTimeFormat("ja-JP", {
  timeZone: "Asia/Tokyo",
  hour: "2-digit",
  minute: "2-digit",
  second: "2-digit",
  hour12: false,
});
function formatJst(julianDate) {
  return jstFormatter.format(Cesium.JulianDate.toDate(julianDate));
}
if (viewer.timeline && typeof viewer.timeline.makeLabel === "function") {
  viewer.timeline.makeLabel = (time) => formatJst(time);
  viewer.timeline.resize();
}
function currentEpoch(julianDate) {
  return Math.floor(Cesium.JulianDate.toDate(julianDate).getTime() / 1000);
}

// -----------------------------------------------------------------------------
// OSM Buildings
// -----------------------------------------------------------------------------
let osmBuildings = null;
try {
  osmBuildings = await Cesium.createOsmBuildingsAsync();
  viewer.scene.primitives.add(osmBuildings);
} catch (error) {
  console.warn("Cesium OSM Buildings could not be loaded:", error);
  showBuildingsEl.checked = false;
  showBuildingsEl.disabled = true;
}

// -----------------------------------------------------------------------------
// Data sources
// -----------------------------------------------------------------------------
const siteDataSource = new Cesium.CustomDataSource("sites");
const routeDataSource = new Cesium.CustomDataSource("routes");
const droneDataSource = new Cesium.CustomDataSource("drones");
const reservationDataSource = new Cesium.CustomDataSource("reservations");
viewer.dataSources.add(siteDataSource);
viewer.dataSources.add(routeDataSource);
viewer.dataSources.add(droneDataSource);
viewer.dataSources.add(reservationDataSource);

// -----------------------------------------------------------------------------
// Sites
// -----------------------------------------------------------------------------
for (const [siteId, site] of Object.entries(dataset.sites)) {
  siteDataSource.entities.add({
    id: `site-${siteId}`,
    name: site.name,
    position: Cesium.Cartesian3.fromDegrees(Number(site.lon), Number(site.lat), 0),
    point: {
      pixelSize: 9,
      color: Cesium.Color.YELLOW,
      outlineColor: Cesium.Color.BLACK,
      outlineWidth: 2,
      heightReference: CLAMP_TO_TERRAIN,
    },
    label: {
      text: site.name,
      font: "14px sans-serif",
      fillColor: Cesium.Color.WHITE,
      outlineColor: Cesium.Color.BLACK,
      outlineWidth: 3,
      style: Cesium.LabelStyle.FILL_AND_OUTLINE,
      pixelOffset: new Cesium.Cartesian2(0, -18),
      heightReference: CLAMP_TO_TERRAIN,
      distanceDisplayCondition: new Cesium.DistanceDisplayCondition(0, 15000),
    },
  });
}

// -----------------------------------------------------------------------------
// Persistent nominal route context
// -----------------------------------------------------------------------------

const routeKeys = new Set();

for (const flight of dataset.flights) {
  const key = `${flight.origin_id}->${flight.destination_id}`;
  if (routeKeys.has(key)) continue;
  routeKeys.add(key);

  const samples = Array.isArray(flight.nominal_trajectory)
    ? flight.nominal_trajectory
    : [];

  const lonLat = [];
  const stride = Math.max(1, Math.floor(samples.length / 160));
  let lastLon = null;
  let lastLat = null;

  for (let i = 0; i < samples.length; i += stride) {
    const s = samples[i];
    const lon = Number(s.lon);
    const lat = Number(s.lat);
    if (!Number.isFinite(lon) || !Number.isFinite(lat)) continue;

    if (
      lastLon !== null &&
      Math.abs(lon - lastLon) < 1e-10 &&
      Math.abs(lat - lastLat) < 1e-10
    ) continue;

    lonLat.push(lon, lat);
    lastLon = lon;
    lastLat = lat;
  }

  if (samples.length) {
    const s = samples[samples.length - 1];
    const lon = Number(s.lon);
    const lat = Number(s.lat);
    if (
      Number.isFinite(lon) &&
      Number.isFinite(lat) &&
      (
        lastLon === null ||
        Math.abs(lon - lastLon) >= 1e-10 ||
        Math.abs(lat - lastLat) >= 1e-10
      )
    ) {
      lonLat.push(lon, lat);
    }
  }

  if (lonLat.length < 4) continue;

  routeDataSource.entities.add({
    id: `route-${flight.origin_id}-${flight.destination_id}`,
    name: `${flight.origin_name} → ${flight.destination_name}`,
    polyline: {
      positions: Cesium.Cartesian3.fromDegreesArray(lonLat),
      width: 3,
      clampToGround: true,
      material: Cesium.Color.YELLOW.withAlpha(0.82),
      distanceDisplayCondition:
        new Cesium.DistanceDisplayCondition(0, 35000),
    },
  });
}

routeDataSource.show = true;

// -----------------------------------------------------------------------------
// Preprocess flight trajectories and dynamic reservation intervals
// -----------------------------------------------------------------------------
function preprocessIntervals(intervals) {
  return intervals.map((iv) => ({
    start: Number(iv.start_epoch_s),
    end: Number(iv.end_epoch_s),
  }));
}
function intervalActive(intervals, epoch) {
  // Intervals are merged and sorted by the generator.
  for (const iv of intervals) {
    if (epoch < iv.start) return false;
    if (epoch < iv.end) return true;
  }
  return false;
}

for (const flight of dataset.flights) {
  flight._actualEpochs = flight.actual_trajectory.map((s) => Number(s.epoch_s));
  for (const cell of flight.reservation.display_cells) {
    cell._intervals = preprocessIntervals(cell.intervals);
  }
}

function sampleAtEpoch(flight, epoch) {
  const epochs = flight._actualEpochs;
  if (!epochs.length || epoch < epochs[0] || epoch > epochs[epochs.length - 1]) {
    return null;
  }
  let lo = 0;
  let hi = epochs.length - 1;
  while (lo <= hi) {
    const mid = (lo + hi) >> 1;
    if (epochs[mid] === epoch) return flight.actual_trajectory[mid];
    if (epochs[mid] < epoch) lo = mid + 1;
    else hi = mid - 1;
  }
  return flight.actual_trajectory[Math.max(0, hi)];
}

// -----------------------------------------------------------------------------
// Index precomputed BST-ID overlap intervals by x/y/h atom.
// -----------------------------------------------------------------------------
function atomKey(x21, y21, h12) {
  return `${x21}/${y21}/${h12}`;
}
const overlapByAtom = new Map();
for (const segment of dataset.overlap_segments ?? []) {
  const cell = segment.cell;
  const key = atomKey(cell.x_index, cell.y_index, cell.h_index);
  if (!overlapByAtom.has(key)) overlapByAtom.set(key, []);
  overlapByAtom.get(key).push({
    start: Number(segment.start_epoch_s),
    end: Number(segment.end_epoch_s),
    coverage: Number(segment.coverage),
    flightIds: segment.flight_ids,
  });
}
for (const segments of overlapByAtom.values()) {
  segments.sort((a, b) => a.start - b.start);
}
// Release the large duplicated top-level array after indexing.
dataset.overlap_segments = null;

function atomCoverageAt(x21, y21, h12, epoch) {
  const segments = overlapByAtom.get(atomKey(x21, y21, h12));
  if (!segments) return 1;
  for (const seg of segments) {
    if (epoch < seg.start) return 1;
    if (epoch < seg.end) return seg.coverage;
  }
  return 1;
}

function cellConflictCoverage(cell, epoch) {
  let max = 1;
  for (const h12 of cell.h12_indices ?? []) {
    max = Math.max(max, atomCoverageAt(cell.x21, cell.y21, h12, epoch));
  }
  return max;
}

function reservationColor(coverage) {
  if (coverage >= 3) return Cesium.Color.RED.withAlpha(0.58);
  if (coverage === 2) return Cesium.Color.ORANGE.withAlpha(0.46);
  return Cesium.Color.CYAN.withAlpha(0.23);
}

// -----------------------------------------------------------------------------
// Moving UAS
// -----------------------------------------------------------------------------
const droneSvg = `
<svg xmlns="http://www.w3.org/2000/svg" width="64" height="64" viewBox="0 0 64 64">
  <g fill="none" stroke="white" stroke-width="5" stroke-linecap="round">
    <path d="M20 20 L44 44 M44 20 L20 44"/>
    <circle cx="14" cy="14" r="8"/><circle cx="50" cy="14" r="8"/>
    <circle cx="14" cy="50" r="8"/><circle cx="50" cy="50" r="8"/>
  </g>
  <circle cx="32" cy="32" r="7" fill="white"/>
</svg>`;
const droneImage = `data:image/svg+xml;charset=UTF-8,${encodeURIComponent(droneSvg)}`;
const flightLabelEntities = [];

for (const flight of dataset.flights) {
  const position = new Cesium.SampledPositionProperty();
  for (const sample of flight.actual_trajectory) {
    position.addSample(
      Cesium.JulianDate.fromIso8601(sample.time_utc),
      Cesium.Cartesian3.fromDegrees(
        Number(sample.lon),
        Number(sample.lat),
        Math.max(0, Number(sample.agl_m))
      )
    );
  }
  position.setInterpolationOptions({
    interpolationDegree: 1,
    interpolationAlgorithm: Cesium.LinearApproximation,
  });

  const actualStart = Cesium.JulianDate.fromIso8601(flight.actual_trajectory[0].time_utc);
  const actualStop = Cesium.JulianDate.fromIso8601(
    flight.actual_trajectory[flight.actual_trajectory.length - 1].time_utc
  );

  const conformanceColor = new Cesium.CallbackProperty((time, result) => {
    const sample = sampleAtEpoch(flight, currentEpoch(time));
    let color = Cesium.Color.WHITE;
    if (sample?.conformance_applicable && sample.inside_reservation === false) {
      color = Cesium.Color.RED;
    } else if (sample && !sample.conformance_applicable) {
      color = Cesium.Color.LIGHTGRAY;
    }
    return color.clone(result);
  }, false);

  const entity = droneDataSource.entities.add({
    id: `drone-${flight.flight_id}`,
    name: `${flight.flight_id}: ${flight.origin_name} → ${flight.destination_name}`,
    availability: new Cesium.TimeIntervalCollection([
      new Cesium.TimeInterval({ start: actualStart, stop: actualStop }),
    ]),
    position,
    billboard: {
      image: droneImage,
      width: 30,
      height: 30,
      color: conformanceColor,
      heightReference: RELATIVE_TO_TERRAIN,
      verticalOrigin: Cesium.VerticalOrigin.CENTER,
      disableDepthTestDistance: 4000,
    },
    label: {
      text: flight.flight_id,
      font: "13px sans-serif",
      fillColor: Cesium.Color.WHITE,
      outlineColor: Cesium.Color.BLACK,
      outlineWidth: 3,
      style: Cesium.LabelStyle.FILL_AND_OUTLINE,
      pixelOffset: new Cesium.Cartesian2(0, -25),
      heightReference: RELATIVE_TO_TERRAIN,
      disableDepthTestDistance: 4000,
    },
    description: new Cesium.CallbackProperty((time) => {
      const sample = sampleAtEpoch(flight, currentEpoch(time));
      if (!sample) return `<b>${flight.flight_id}</b>`;
      let state = "N/A (takeoff/landing)";
      if (sample.conformance_applicable) {
        state = sample.inside_reservation ? "INSIDE" : `OUTSIDE (${sample.violation})`;
      }
      return `
        <table>
          <tr><td><b>Flight</b></td><td>${flight.flight_id}</td></tr>
          <tr><td><b>Route</b></td><td>${flight.origin_name} → ${flight.destination_name}</td></tr>
          <tr><td><b>Target AGL</b></td><td>${Number(flight.target_agl_m).toFixed(0)} m</td></tr>
          <tr><td><b>Phase</b></td><td>${sample.phase}</td></tr>
          <tr><td><b>AGL</b></td><td>${Number(sample.agl_m).toFixed(1)} m</td></tr>
          <tr><td><b>Speed</b></td><td>${Number(sample.horizontal_speed_mps).toFixed(1)} m/s</td></tr>
          <tr><td><b>Heading error</b></td><td>${Number(sample.heading_error_deg).toFixed(1)}°</td></tr>
          <tr><td><b>Cross-track</b></td><td>${Number(sample.cross_track_m).toFixed(1)} m</td></tr>
          <tr><td><b>Conformance</b></td><td><b>${state}</b></td></tr>
        </table>`;
    }, false),
  });
  flightLabelEntities.push(entity);
}

// -----------------------------------------------------------------------------
// Dynamic reservation layer.
// The generator has already compressed 1-Hz sample reservations into merged
// validity intervals for each z21 display cell.
// -----------------------------------------------------------------------------
const RENDER_STEP_S = 5;
let lastRenderBucket = null;
let currentReservationStats = {
  activeReservations: 0,
  renderedCells: 0,
  overlapCells: 0,
  maxCoverage: 0,
};

function rebuildReservationLayer(epoch, force = false) {
  const bucket = Math.floor(epoch / RENDER_STEP_S);
  const modeSignature = `${showReservationsEl.checked}/${overlapOnlyEl.checked}`;
  const signature = `${bucket}/${modeSignature}`;
  if (!force && signature === lastRenderBucket) return;
  lastRenderBucket = signature;

  reservationDataSource.entities.suspendEvents();
  reservationDataSource.entities.removeAll();

  if (!showReservationsEl.checked) {
    reservationDataSource.entities.resumeEvents();
    currentReservationStats = {
      activeReservations: 0,
      renderedCells: 0,
      overlapCells: 0,
      maxCoverage: 0,
    };
    return;
  }

  // Aggregate coincident cells at the same nominal AGL band so that two
  // reservations do not simply paint the same rectangle twice.
  const groups = new Map();
  const activeFlights = new Set();

  for (const flight of dataset.flights) {
    const bottom = Number(flight.reservation.agl_bottom_m);
    const top = Number(flight.reservation.agl_top_m);

    for (const cell of flight.reservation.display_cells) {
      if (!intervalActive(cell._intervals, epoch)) continue;
      activeFlights.add(flight.flight_id);

      const key = `${cell.x21}/${cell.y21}/${bottom}/${top}`;
      let rec = groups.get(key);
      if (!rec) {
        rec = {
          cell,
          bottom,
          top,
          flightIds: new Set(),
          hIndices: new Set(),
        };
        groups.set(key, rec);
      }
      rec.flightIds.add(flight.flight_id);
      for (const h of cell.h12_indices ?? []) rec.hIndices.add(Number(h));
    }
  }

  let renderedCells = 0;
  let overlapCells = 0;
  let maxCoverage = 0;

  for (const rec of groups.values()) {
    const cell = rec.cell;
    let coverage = 1;
    for (const h12 of rec.hIndices) {
      coverage = Math.max(
        coverage,
        atomCoverageAt(Number(cell.x21), Number(cell.y21), h12, epoch)
      );
    }
    if (coverage >= 2) overlapCells += 1;
    maxCoverage = Math.max(maxCoverage, coverage);
    if (overlapOnlyEl.checked && coverage < 2) continue;

    const ids = Array.from(rec.flightIds).sort();
    reservationDataSource.entities.add({
      id: `res-${cell.x21}-${cell.y21}-${rec.bottom}-${rec.top}`,
      name: `Reservation ${rec.bottom.toFixed(0)}-${rec.top.toFixed(0)} m AGL`,
      rectangle: {
        coordinates: Cesium.Rectangle.fromDegrees(
          Number(cell.bounds.west),
          Number(cell.bounds.south),
          Number(cell.bounds.east),
          Number(cell.bounds.north)
        ),
        // Exact intended operational band: 30 m route => 20-40 m AGL, etc.
        height: rec.bottom,
        heightReference: RELATIVE_TO_TERRAIN,
        extrudedHeight: rec.top,
        extrudedHeightReference: RELATIVE_TO_TERRAIN,
        material: reservationColor(coverage),
        outline: coverage >= 2,
        outlineColor: coverage >= 3 ? Cesium.Color.WHITE : Cesium.Color.ORANGE,
      },
      description: `
        <table>
          <tr><td><b>Flights</b></td><td>${ids.join(", ")}</td></tr>
          <tr><td><b>Coverage</b></td><td>${coverage}</td></tr>
          <tr><td><b>XY atom</b></td><td>${cell.x21}/${cell.y21} (z${dataset.metadata.reservation.margin_xy_zoom})</td></tr>
          <tr><td><b>Reserved AGL</b></td><td>${rec.bottom.toFixed(1)}-${rec.top.toFixed(1)} m</td></tr>
          <tr><td><b>H indices</b></td><td>${Array.from(rec.hIndices).sort((a,b)=>a-b).join(", ")}</td></tr>
          <tr><td><b>Time rule</b></td><td>expected passage ±5 min</td></tr>
        </table>`,
    });
    renderedCells += 1;
  }

  reservationDataSource.entities.resumeEvents();
  currentReservationStats = {
    activeReservations: activeFlights.size,
    renderedCells,
    overlapCells,
    maxCoverage,
  };
}

// -----------------------------------------------------------------------------
// Status
// -----------------------------------------------------------------------------
let lastStatusEpoch = null;
function updateStatus(epoch) {
  if (epoch === lastStatusEpoch) return;
  lastStatusEpoch = epoch;

  let activeFlights = 0;
  let applicableFlights = 0;
  let insideFlights = 0;
  let outsideFlights = 0;
  const outsideIds = [];

  for (const flight of dataset.flights) {
    const sample = sampleAtEpoch(flight, epoch);
    if (!sample) continue;
    activeFlights += 1;
    if (!sample.conformance_applicable) continue;
    applicableFlights += 1;
    if (sample.inside_reservation) insideFlights += 1;
    else {
      outsideFlights += 1;
      outsideIds.push(flight.flight_id);
    }
  }

  const outsideText = outsideIds.length
    ? `<br><span style="color:#ff7777">Outside: ${outsideIds.join(", ")}</span>`
    : "";

  scenarioStatusEl.innerHTML = `
    <b>JST ${formatJst(viewer.clock.currentTime)}</b><br>
    Flying: ${activeFlights}
    &nbsp; | &nbsp; Cruise monitored: ${applicableFlights}
    &nbsp; | &nbsp; Inside: ${insideFlights}
    &nbsp; | &nbsp; Outside: ${outsideFlights}<br>
    Active route reservations: ${currentReservationStats.activeReservations}
    &nbsp; | &nbsp; Rendered cells: ${currentReservationStats.renderedCells}<br>
    Overlap cells: ${currentReservationStats.overlapCells}
    &nbsp; | &nbsp; Max coverage: ${currentReservationStats.maxCoverage}
    ${outsideText}`;
}

// -----------------------------------------------------------------------------
// Camera and controls
// -----------------------------------------------------------------------------
function setAizuView() {
  viewer.camera.setView({
    destination: Cesium.Cartesian3.fromDegrees(139.9405, 37.5060, 5200),
    orientation: {
      heading: Cesium.Math.toRadians(8),
      pitch: Cesium.Math.toRadians(-55),
      roll: 0,
    },
  });
}
setAizuView();

speedSelectEl.addEventListener("change", () => {
  viewer.clock.multiplier = Number(speedSelectEl.value);
});
showReservationsEl.addEventListener("change", () => {
  rebuildReservationLayer(currentEpoch(viewer.clock.currentTime), true);
});
overlapOnlyEl.addEventListener("change", () => {
  rebuildReservationLayer(currentEpoch(viewer.clock.currentTime), true);
});
showLabelsEl.addEventListener("change", () => {
  for (const entity of flightLabelEntities) {
    if (entity.label) entity.label.show = showLabelsEl.checked;
  }
});
showRoutesEl.addEventListener("change", () => {
  routeDataSource.show = showRoutesEl.checked;
});

baseMapSelectEl.addEventListener("change", () => {
  gsiBaseLayer = addGsiBaseMap(baseMapSelectEl.value);
});

showBuildingsEl.addEventListener("change", () => {
  if (osmBuildings) osmBuildings.show = showBuildingsEl.checked;
});

document.getElementById("homeView").addEventListener("click", setAizuView);

viewer.clock.onTick.addEventListener((clock) => {
  const epoch = currentEpoch(clock.currentTime);
  rebuildReservationLayer(epoch);
  updateStatus(epoch);
});

const initialEpoch = currentEpoch(viewer.clock.currentTime);
rebuildReservationLayer(initialEpoch, true);
updateStatus(initialEpoch);
console.log("BST-ID dynamic UTM scenario loaded", dataset.summary);


// -----------------------------------------------------------------------------
// Optional static Tsuruga Castle layer
// -----------------------------------------------------------------------------

const castleToolbar = document.getElementById("castleToolbar");

const tsurugaDataUrl =
  window.TSURUGA_BSTID_DATA_URL ?? "../tsuruga_native_1m_cesium.json";

try {
  const { installTsurugaBstidLayer } =
    await import("./tsuruga_layer.js");

  await installTsurugaBstidLayer(viewer, {
    dataUrl: tsurugaDataUrl,
    toolbar: castleToolbar,
  });
} catch (error) {
  console.warn("Tsuruga Castle layer unavailable:", error);

  if (castleToolbar) {
    castleToolbar.innerHTML = `
      <div class="title">Tsuruga Castle BST-ID Object</div>
      <div class="status-box error-box">
        Castle layer could not be loaded.<br>
        ${tsurugaDataUrl}
      </div>
    `;
  }
}

