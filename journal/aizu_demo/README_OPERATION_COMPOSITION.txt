BST-ID UTM operation-composition viewer

Replace:
  journal/aizu_demo/generate_utm_scenario.py
  journal/aizu_demo/viewer/app.js

Then regenerate:
  cd C:\work\bst-id
  python journal\aizu_demo\generate_utm_scenario.py

Serve:
  cd journal\aizu_demo
  py -m http.server 8000

Open:
  http://localhost:8000/viewer/

The viewer executes at each Cesium time:
  TimeSlice(t)
    -> Union(active XYH reservation atoms)
    -> Normalize(XYH)
    -> Project H to the operational AGL band
    -> Normalize(XY for display)

Sibling prefixes are merged only when both the operational altitude band
and the exact active-flight set are identical. Therefore cyan reservation
regions do not merge across orange/red overlap boundaries.

Use the 'TimeSlice -> Normalize' checkbox to compare normalized display
against the unmerged z21 projection.
