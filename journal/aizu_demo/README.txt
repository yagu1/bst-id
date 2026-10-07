Multi-speed update

Cruise speed per flight:
  5 m/s, 8 m/s, or 10 m/s

The 34 flights are assigned approximately evenly and shuffled
deterministically using the scenario random seed.

Speed-dependent behavior:
  5 m/s  calmer disturbance, 70 m look-ahead, approx +/-6% speed variation
  8 m/s  medium disturbance, 112 m look-ahead, approx +/-8% speed variation
 10 m/s  current-like disturbance, 140 m look-ahead, approx +/-10% speed variation

Climb/descent remain 3 m/s.

Replace:
  journal/aizu_demo/generate_utm_scenario.py
  journal/aizu_demo/viewer/app.js

Then regenerate utm_scenario.json and reload the browser with Ctrl+F5.
