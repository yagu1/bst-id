# BST-ID whole-scenario compression analysis

## Evaluation scope

- Flights: 34
- Fixed working grid: X z21, Y z21, H z12
- Semantic attribute: exact active reservation-flight set
- Merge rule: prefixes are mergeable only when their semantic attribute is identical.
- Attribute comparison model: 34-bit flight-set mask per record.

The main comparison isolates the effect of BST-ID prefix normalization:

1. **Fixed grid + dyadic T**: XYH remains fixed at z21/z21/H12; exact time intervals are represented by BST-ID-compatible dyadic T prefixes.
2. **BST-ID normalized XYHT**: the same semantic data are normalized by independent binary sibling merging along X, Y, H, and T.

The arbitrary exact-interval baseline is also reported as a strong practical reference. It can encode any start/end pair directly, whereas BST-ID T prefixes are dyadic.

## Main result

- Fixed-grid dyadic records: **932,090**
- Normalized BST-ID records: **240,360**
- Record reduction: **74.21%**
- Record compression factor: **3.878x**
- Key-bit reduction: **74.47%**
- Key-bit compression factor: **3.916x**
- Total-bit reduction including the controlled attribute mask: **74.40%**

## Why records were reduced

Binary sibling merges:

- X-axis: **332,006**
- Y-axis: **227,458**
- H-axis: **132,266**
- T-axis: **0**

Final normalized records with at least one spatial axis coarsened:
**216,818/240,360 (90.21%)**

Final records with different X/Y zooms (direct evidence of rectangular, non-square horizontal prefixes):
**181,993/240,360 (75.72%)**

## Output files

- `compression_summary.json`: complete machine-readable summary
- `compression_summary.csv`: representation-level comparison
- `axis_merge_summary.csv`: merge counts by X/Y/H/T
- `zoom_distribution.csv`: final zoom-vector distribution
- `attribute_cardinality_distribution.csv`: ordinary vs overlap-region statistics
- `flight_summary.csv`: per-flight context including 5/8/10 m/s cruise speed
- `fixed_interval_records.csv.gz`: full fixed-grid exact-interval records
- `normalized_xyht_records.csv.gz`: full normalized XYHT prefix records

## Interpretation note

This analysis does **not** claim a globally minimum rectangle cover. It applies a deterministic X->Y->H->T binary sibling-normalization policy while preserving the exact semantic attribute. The comparison therefore evaluates the actual hierarchical representation produced by the stated policy, not an unconstrained optimal packing.
