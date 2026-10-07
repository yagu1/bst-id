// tsuruga_layer.js
//
// Drop-in CesiumJS module for the existing BST-ID UTM viewer.
//
// Usage from the existing UTM app.js, AFTER viewer and toolbar exist:
//
//   const { installTsurugaBstidLayer } = await import("./tsuruga_layer.js");
//   await installTsurugaBstidLayer(viewer, {
//     dataUrl: "../tsuruga_native_cesium.json",
//     toolbar,
//   });
//
// Tsuruga Castle uses one common local-U vertical datum.
// We therefore sample terrain height ONCE at the castle anchor and render all
// BST-ID cells at absolute ellipsoid heights:
//
//   display height = anchor terrain height + local U
//
// This intentionally avoids RELATIVE_TO_TERRAIN for each rectangle, because
// large normalized cells can span terrain-height variation and otherwise
// acquire an apparent vertical distortion that is not present in BST-ID.

function cellDescription(prefix, label) {
  return `
    <table>
      <tr><td><b>Layer</b></td><td>${label}</td></tr>
      <tr><td><b>X</b></td><td>z${prefix.zx} / ${prefix.ix}</td></tr>
      <tr><td><b>Y</b></td><td>z${prefix.zy} / ${prefix.iy}</td></tr>
      <tr><td><b>H</b></td><td>z${prefix.zh} / ${prefix.ih}</td></tr>
      <tr><td><b>Local U</b></td>
          <td>${prefix.bottom_u.toFixed(2)}–${prefix.top_u.toFixed(2)} m</td></tr>
    </table>
  `;
}

function addPrefixEntity(ds, prefix, options) {
  const bottom = Math.max(0.0, Number(prefix.bottom_u));
  const top = Math.max(bottom + 0.05, Number(prefix.top_u));
  const baseHeight = Number(options.baseHeight ?? 0.0);

  const absoluteBottom = baseHeight + bottom;
  const absoluteTop = baseHeight + top;

  return ds.entities.add({
    id: options.id,
    name: options.name,
    rectangle: {
      coordinates: Cesium.Rectangle.fromDegrees(
        prefix.west,
        prefix.south,
        prefix.east,
        prefix.north
      ),
      height: absoluteBottom,
      heightReference: Cesium.HeightReference.NONE,
      extrudedHeight: absoluteTop,
      extrudedHeightReference: Cesium.HeightReference.NONE,
      material: options.color,
      outline: Boolean(options.outline),
      outlineColor: options.outlineColor ?? Cesium.Color.WHITE,
    },
    description: cellDescription(prefix, options.label),
  });
}

async function resolveAnchorTerrainHeight(viewer, dataset, explicitBaseHeight) {
  if (Number.isFinite(Number(explicitBaseHeight))) {
    return Number(explicitBaseHeight);
  }

  const lon = Number(dataset.metadata.anchor_longitude_deg);
  const lat = Number(dataset.metadata.anchor_latitude_deg);
  const anchor = Cesium.Cartographic.fromDegrees(lon, lat, 0.0);

  // Preferred path: sample the terrain provider directly.  This deliberately
  // ignores 3D Tiles/buildings and obtains the terrain surface only.
  try {
    if (viewer.terrainProvider && Cesium.sampleTerrainMostDetailed) {
      const sampled = await Cesium.sampleTerrainMostDetailed(
        viewer.terrainProvider,
        [anchor]
      );
      const h = sampled?.[0]?.height;
      if (Number.isFinite(h)) {
        return h;
      }
    }
  } catch (err) {
    console.warn(
      "Tsuruga BST-ID: terrain sampling at anchor failed; trying globe cache.",
      err
    );
  }

  // Fallback: use any terrain height already available in the globe cache.
  try {
    const h = viewer.scene?.globe?.getHeight(anchor);
    if (Number.isFinite(h)) {
      return h;
    }
  } catch (err) {
    console.warn("Tsuruga BST-ID: globe.getHeight fallback failed.", err);
  }

  console.warn(
    "Tsuruga BST-ID: no terrain height available at anchor; " +
      "using ellipsoid height 0 m as the local-U base."
  );
  return 0.0;
}

function setSourceVisible(ds, visible) {
  ds.show = Boolean(visible);
}

export async function installTsurugaBstidLayer(
  viewer,
  {
    dataUrl = "../tsuruga_native_cesium.json",
    toolbar = document.getElementById("toolbar"),
    // Optional override in metres above the WGS84 ellipsoid.
    // Leave undefined to sample terrain once at the castle anchor.
    anchorBaseHeightM = undefined,
  } = {}
) {
  let dataset;
  const response = await fetch(dataUrl);
  if (!response.ok) {
    throw new Error(
      `Failed to load Tsuruga BST-ID data: HTTP ${response.status}`
    );
  }
  dataset = await response.json();

  const baseHeight = await resolveAnchorTerrainHeight(
    viewer,
    dataset,
    anchorBaseHeightM
  );

  console.log(
    `Tsuruga BST-ID vertical datum: anchor terrain ellipsoid height = ` +
      `${baseHeight.toFixed(3)} m`
  );

  const bitmapDS = new Cesium.CustomDataSource("tsuruga-bitmap");
  const normalizedDS = new Cesium.CustomDataSource("tsuruga-normalized");
  const specialDS = new Cesium.CustomDataSource("tsuruga-special");

  viewer.dataSources.add(bitmapDS);
  viewer.dataSources.add(normalizedDS);
  viewer.dataSources.add(specialDS);

  // Fixed-grid fine BST-ID occupancy.
  bitmapDS.entities.suspendEvents();
  dataset.bitmap.forEach((cell, i) => {
    addPrefixEntity(bitmapDS, cell, {
      id: `tsuruga-bitmap-${i}`,
      name: "Fine BST-ID occupancy cell",
      label: "Fine BST-ID / bitmap-equivalent",
      baseHeight,
      color: Cesium.Color.LIGHTGRAY.withAlpha(0.30),
      outline: false,
    });
  });
  bitmapDS.entities.resumeEvents();

  // Normalized hierarchical BST-ID.
  normalizedDS.entities.suspendEvents();
  dataset.normalized.forEach((cell, i) => {
    const coarsening =
      (dataset.metadata.fine_zoom.x - cell.zx) +
      (dataset.metadata.fine_zoom.y - cell.zy) +
      (dataset.metadata.fine_zoom.h - cell.zh);

    const alpha = 0.50;

    addPrefixEntity(normalizedDS, cell, {
      id: `tsuruga-normalized-${i}`,
      name: "Normalized BST-ID prefix",
      label: "Normalized hierarchical BST-ID",
      baseHeight,
      color: Cesium.Color.DODGERBLUE.withAlpha(alpha),
      outline: coarsening >= 4,
      outlineColor: Cesium.Color.CYAN,
    });
  });
  normalizedDS.entities.resumeEvents();

  // Maximum inscribed cell.
  const inner = dataset.maximum_inscribed;
  const innerEntity = addPrefixEntity(specialDS, inner, {
    id: "tsuruga-maximum-inscribed",
    name: "Maximum inscribed BST-ID cell",
    label: "Maximum inscribed BST-ID cell",
    baseHeight,
    color: Cesium.Color.LIME.withAlpha(0.58),
    outline: true,
    outlineColor: Cesium.Color.LIME,
  });

  // Minimum enclosing cell / cover.
  const outerEntities = [];
  dataset.minimum_enclosing.cells.forEach((cell, i) => {
    outerEntities.push(
      addPrefixEntity(specialDS, cell, {
        id: `tsuruga-minimum-enclosing-${i}`,
        name:
          dataset.minimum_enclosing.mode === "single_cell"
            ? "Minimum enclosing BST-ID cell"
            : "Minimum enclosing BST-ID cover",
        label:
          dataset.minimum_enclosing.mode === "single_cell"
            ? "Minimum enclosing BST-ID cell"
            : "Minimum valid enclosing BST-ID cover",
        baseHeight,
        color: Cesium.Color.RED.withAlpha(0.08),
        outline: true,
        outlineColor: Cesium.Color.RED,
      })
    );
  });

  // Initial state: normalized hierarchy only.
  bitmapDS.show = false;
  normalizedDS.show = true;
  innerEntity.show = true;
  outerEntities.forEach((e) => (e.show = true));

  // Add a compact section to the existing UTM toolbar.
  if (toolbar) {
    const section = document.createElement("div");
    section.id = "tsurugaBstidControls";
    section.style.marginTop = "10px";
    section.style.paddingTop = "8px";
    section.style.borderTop = "1px solid rgba(255,255,255,0.28)";

    section.innerHTML = `
      <div style="font-weight:600;margin-bottom:5px;">Tsuruga Castle BST-ID</div>
      <div style="display:flex;flex-wrap:wrap;gap:7px 12px;align-items:center;">
        <label>
          <input id="tsurugaBitmap" type="checkbox">
          Fine bitmap
        </label>
        <label>
          <input id="tsurugaNormalized" type="checkbox" checked>
          Normalize
        </label>
        <label>
          <input id="tsurugaInner" type="checkbox" checked>
          Max inscribed
        </label>
        <label>
          <input id="tsurugaOuter" type="checkbox" checked>
          Min enclosing
        </label>
        <button id="tsurugaView" type="button">Castle view</button>
      </div>
      <div id="tsurugaStatus"
           style="margin-top:5px;font-size:12px;line-height:1.35;">
        Fine ${dataset.metadata.fine_cells.toLocaleString()} cells →
        ${dataset.metadata.normalized_prefixes.toLocaleString()} prefixes
        (${dataset.metadata.reduction_percent.toFixed(1)}% reduction)<br>
        Fixed vertical datum: anchor terrain ${baseHeight.toFixed(2)} m
        (WGS84 ellipsoid)
      </div>
    `;

    toolbar.appendChild(section);

    document
      .getElementById("tsurugaBitmap")
      .addEventListener("change", (e) => setSourceVisible(bitmapDS, e.target.checked));

    document
      .getElementById("tsurugaNormalized")
      .addEventListener("change", (e) => setSourceVisible(normalizedDS, e.target.checked));

    document
      .getElementById("tsurugaInner")
      .addEventListener("change", (e) => {
        innerEntity.show = e.target.checked;
      });

    document
      .getElementById("tsurugaOuter")
      .addEventListener("change", (e) => {
        outerEntities.forEach((entity) => {
          entity.show = e.target.checked;
        });
      });

    document
      .getElementById("tsurugaView")
      .addEventListener("click", () => {
        flyToCastle(viewer, dataset, baseHeight);
      });
  }

  flyToCastle(viewer, dataset, baseHeight);

  console.log(
    `Tsuruga BST-ID loaded: fine=${dataset.metadata.fine_cells}, ` +
    `normalized=${dataset.metadata.normalized_prefixes}, ` +
    `outer=${dataset.minimum_enclosing.cells.length}`
  );

  return {
    dataset,
    bitmapDS,
    normalizedDS,
    specialDS,
    baseHeight,
  };
}

function flyToCastle(viewer, dataset, baseHeight = 0.0) {
  const lon = Number(dataset.metadata.anchor_longitude_deg);
  const lat = Number(dataset.metadata.anchor_latitude_deg);

  viewer.camera.flyTo({
    destination: Cesium.Cartesian3.fromDegrees(lon, lat, baseHeight + 180.0),
    orientation: {
      heading: Cesium.Math.toRadians(20.0),
      pitch: Cesium.Math.toRadians(-48.0),
      roll: 0.0,
    },
    duration: 1.4,
  });
}
