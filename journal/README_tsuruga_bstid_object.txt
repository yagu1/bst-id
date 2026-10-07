Tsuruga Castle point cloud -> filled BST-ID object
================================================

Files
-----
tsuruga_ply_to_bstid_object.py

Purpose
-------
Takes the already aligned Tsuruga Castle PLY:

  +X = East
  +Y = North
  +Z = Up
  anchor = selected castle/ground boundary point

and creates:

1. surface voxels;
2. a filled 3-D object;
3. fine XYH BST-ID cells;
4. a normalized BST-ID region.

Install
-------
In the BST-ID virtual environment:

  pip install open3d scipy pyproj numpy

Run from the BST-ID repository root so that:

  from bst_id.region_algebra import ...

works.

Recommended first run
---------------------

  python tsuruga_ply_to_bstid_object.py ^
      tsuruga_aligned.ply ^
      --transform-json tsuruga_enu_anchor_transform.json ^
      --voxel-size 1.0 ^
      --fill-mode directional ^
      --close-iterations 1 ^
      --zx 25 --zy 25 --zh 15 ^
      --verify ^
      --out-prefix tsuruga_1m

PowerShell uses ^ only in cmd.exe. In PowerShell either use one line or use
the backtick (`) for line continuation.

Fill modes
----------
directional
    Implements the proposed geometric rule:
    an empty voxel is filled when shell voxels exist to its west, east,
    and above.

    Add --require-ns to additionally require north/south enclosure.

flood
    Uses a 3-D exterior flood-fill interpretation:
    empty cells not connected to the outside are filled.

Suggested workflow
------------------
First compare:

  --fill-mode directional

and:

  --fill-mode flood

by opening:

  *_solid_voxels.ply

in Open3D / CloudCompare.

Use the one whose object interpretation is appropriate and document that
rule explicitly in the manuscript.

Outputs
-------
<out>_shell_voxels.ply
<out>_solid_voxels.ply
<out>_normalized_prefixes.csv
<out>_stats.json

The JSON includes:
- shell / interior / solid voxel counts;
- fine BST-ID cell count;
- normalized prefix count;
- reduction percentage;
- per-axis coarsening counts;
- anisotropic-prefix count;
- dominant (zx,zy,zh) vectors;
- optional exact fine-grid verification.

Altitude
--------
The current transform JSON has no absolute altitude. Therefore, unless
--anchor-alt is supplied, the script uses:

  h = local U, with anchor h = 0 m.

This is suitable for the first local-object experiment. If the final paper
requires globally meaningful BST-ID altitude, rerun with the physical anchor
altitude:

  --anchor-alt <metres>

Normalization
-------------
The default is:

  --normalizer indexed-batch

because a castle solid may contain hundreds of thousands of fine cells.
The transparent reference normalize() has a high conservative worst-case cost
and should only be used on small test cases:

  --normalizer reference

Exactness check
---------------
--verify expands the normalized region back to the fine XYH working zoom and
checks exact equality with the fine BST-ID occupancy.

This can consume substantial memory for a large object but is desirable for
the final reported experiment.
