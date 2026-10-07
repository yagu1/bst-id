Tsuruga Castle BST-ID visualizer
================================

Purpose
-------
Displays four related representations in one Open3D visualizer:

1. Bitmap/fixed-grid occupancy
2. Normalized hierarchical BST-ID fill
3. Maximum inscribed BST-ID cell
4. Minimum enclosing BST-ID cell

The maximum inscribed cell is computed from the original fine BST-ID
occupancy with a 3-D summed-volume table.  It is not merely the largest
prefix already present in the normalized output.

The minimum enclosing cell is the component-wise longest-common-prefix
(LCP) cell containing all fine occupied cells.

Install
-------
  pip install open3d pyproj numpy

Run
---
Place the script in the BST-ID repository/work directory together with:

  tsuruga_1m_solid_voxels.ply
  tsuruga_1m_normalized_prefixes.csv
  tsuruga_1m_stats.json

Then:

  python visualize_tsuruga_bstid.py

Viewer keys
-----------
  B : bitmap / fixed-grid representation
  H : normalized hierarchical BST-ID representation
  I : toggle maximum inscribed BST-ID cell
  E : toggle minimum enclosing BST-ID cell
  A : show all
  C : hide bitmap / hierarchy
  P : save screenshot
  Q : close

Interpretation
--------------
Maximum inscribed BST-ID cell:
  largest single BST-ID cell, measured by fine-atom capacity, whose entire
  fine descendant set is occupied by the object.

Minimum enclosing BST-ID cell:
  smallest single anisotropic BST-ID cell whose x, y, and h prefixes contain
  all fine occupied cells.

Thus:

  C_in  subset_of  R_object  subset_of  C_out

Outputs
-------
  tsuruga_bstid_visualizer_report.json
  tsuruga_hierarchical_bstid.obj
  tsuruga_bstid_visualizer.png  (after pressing P)

The hierarchy is also exported as a triangle-mesh OBJ because Windows 3D
viewers often handle OBJ meshes more reliably than point-cloud-only PLY files.
