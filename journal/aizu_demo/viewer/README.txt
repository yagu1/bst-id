BST-ID UTM + Tsuruga Castle integrated viewer v4
===================================================

This version addresses the black/no-map display directly.

Cause
-----
The actual old UTM app relied on Cesium's default imagery provider.  In the
current integrated environment, the terrain and Cesium entities can still
render even when that default imagery fails or is unavailable.  The result is
a dark/black globe with drones and reservation cells floating on it.

Fix
---
v4 explicitly loads GSI tiles and no longer depends on the Cesium default
basemap.

The UTM panel now includes:

  Base map
    GSI standard
    GSI aerial

The viewer is created with baseLayerPicker=false and GSI imagery is added with
Cesium.UrlTemplateImageryProvider.

The UTM and Castle controls are also separated into TWO visible panels:

  top-left  : UTM
  top-right : Tsuruga Castle BST-ID Object

so the castle controls are no longer hidden below a scroll area.

Files
-----
index.html
app.js
tsuruga_layer.js
config.local.example.js
README.txt

Installation
------------
Replace/add in:

  journal/aizu_demo/viewer/

  index.html
  app.js
  tsuruga_layer.js

Keep your own:
  config.local.js

It should contain:

  window.CESIUM_ION_TOKEN = "...";
  window.TSURUGA_BSTID_DATA_URL =
    "../tsuruga_native_1m_cesium.json";

Expected data one level above viewer/
-------------------------------------
utm_scenario.json
tsuruga_native_1m_cesium.json

Run
---
  cd C:\work\bst-id\journal\aizu_demo
  py -m http.server 8000

Then hard-refresh the browser:

  Ctrl + F5

Open:
  http://localhost:8000/viewer/

Expected
--------
1. GSI standard map should be visible immediately.
2. UTM controls appear at top-left.
3. Castle controls appear at top-right.
4. Aizu view shows persistent yellow nominal route lines plus dynamic UAS /
   reservation layers.
5. Castle view frames the normalized castle data using viewer.flyTo().
6. Min enclosing is off by default.
7. Fine bitmap is lazy-built.

Important
---------
The Python HTTP server log will NOT show GSI tile requests because those are
made by the browser directly to cyberjapandata.gsi.go.jp.  Use the browser
DevTools Network/Console if GSI tiles still fail.

If the map is still black after Ctrl+F5:
- open browser DevTools -> Console;
- check for errors mentioning cyberjapandata.gsi.go.jp, imagery, CORS, or
  mixed-content policy;
- switch Base map from GSI standard to GSI aerial and back once.
