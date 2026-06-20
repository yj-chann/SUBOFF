# SUBOFF Bow H-O Grid Points Generator

This project contains MATLAB scripts for generating the point coordinates of a
three-dimensional H-O type structured grid around the bow of the DARPA SUBOFF
axisymmetric hull, or another body of revolution if the geometry helper
functions are replaced.

The main output is an OpenFOAM v2312 `constant/polyMesh/points`-format file.
The code also writes several Tecplot `.plt` files that are useful for checking
the 2D section, the H-type RECT surface patch, and the O-type CY surface before generating
the points in OpenFOAM.

## Grid Idea

The grid follows the topology shown in the supplied reference figure (Waiting for completed ...):

- An inner H-type block is placed at the bow tip. This avoids the polar
  centerline singularity that would occur if a pure O-grid were collapsed at
  the nose.
- The H-block is smoothly blended into an O-type surface grid around the
  circular body.
- Points are then marched in the outward normal direction from the body surface
  to form the three-dimensional volume coordinates.

In the code this is done in two parts:

- `RECT_*` arrays define the rectangular/H-type bow cap used by the H-block.
- `CY_*` arrays define the cylindrical/O-type surface around the hull.

The final OpenFOAM point list is written in the same ordering expected by the
matching OpenFOAM topology files. This script generates only `points`; the
corresponding `faces`, `owner`, `neighbour`, and `boundary` files must use the
same indexing convention.

## Repository Layout

```text
02_grid_use/
├── README.md
└── file0_code/
    ├── grid.m                  # Main script: builds surface grids and writes points
    ├── splitEdge.m             # OpenFOAM-like multi-grading edge splitter
    ├── Get_R.m                 # SUBOFF radius R(x)
    ├── Get_K.m                 # SUBOFF slope dR/dx
    ├── Get_S.m                 # Hull arc length S(x)
    ├── Get_x_from_S.m          # Inverse mapping from arc length to x
    ├── Get_x_from_R.m          # Inverse mapping from radius to x, bow only
    ├── Get_P.m                 # Normal-offset point from hull surface
    └── Tecplot_InputFiles/     # Tecplot outputs written by grid.m
```

## Quick Start

1. Open MATLAB.
2. Change the MATLAB working directory to `file0_code`.
3. Make sure `Tecplot_InputFiles` exists.
4. Run:

   ```matlab
   grid
   ```

5. Inspect the Tecplot files first, especially `suboff_mesh_2d.plt`,
   `RECT_surface.plt`, and `CY_surface.plt`.
6. Move or copy the generated `points` file into the OpenFOAM case path:

   ```text
   09_3d_startup/constant/polyMesh/points
   ```

The script opens output files using relative paths, so the run directory
matters. If you run `grid.m` from another folder, MATLAB will write `points` and
`Tecplot_InputFiles/*` relative to that current folder.

## Important Resolution Variables

The main resolution controls are at the top of `grid.m`:

```matlab
NJ1=150; NJ2=450; NJ3=300; NK=640;
NL=NK/4; NJ=NJ1+NJ2+NJ3;
NI=400; NI_ADD=40; NJ_ADD=1;
```

Their meanings are:

| Variable | Meaning |
| --- | --- |
| `NK` | Number of circumferential points around the O-grid surface, excluding the duplicated closing point in Tecplot output. Must be divisible by 8 because the H-O transition uses quarter and eighth sectors. |
| `NL` | Derived square/H-block side control, `NL = NK/4`. Do not set independently unless the indexing logic is rewritten. |
| `NJ1` | Number of surface cells in the transition from square H-block to circular O-grid. |
| `NJ2` | Number of surface cells in the next streamwise hull section. |
| `NJ3` | Number of surface cells in the final streamwise hull section. |
| `NJ` | Total streamwise surface cells, `NJ1 + NJ2 + NJ3`. |
| `NI` | Number of wall-normal cells actually written into the OpenFOAM `points` file. |
| `NI_ADD` | Extra wall-normal cells used for `suboff_mesh_2d_add.plt` only. They are not written to the OpenFOAM `points` file, but used in `05_ini_2d_RANS` to calculate flow in an extended domain in order to obtain boundary value for the smaller 3D domain. |
| `NJ_ADD` | Extra streamwise cell used for `suboff_mesh_2d_add.plt` only, similar to `NI_ADD`. |

With the default values, the point count is:

```text
P_NUM = (NI+1)*(NJ+1)*NK + (NL-1)^2*(NI+1)
      = 241,370,321 points
```

This is a very large ASCII `points` file. Use a much smaller trial resolution
before running a production case, for example by reducing `NK`, `NI`, and the
`NJ*` values while preserving the required divisibility and indexing
relationships.

## Wall-Normal Spacing

The wall-normal coordinate is built with:

```matlab
YL = 0.55;
YlenRatios  = [0.005 , 0.18 , 1 , 5.75, 10];
YcellRatios = [10 , 70 , 70 , 50 , 20];
YexpRatios  = [1 , 10.7, 3.9, 10.9 , 1];
Y = splitEdge(YL, NN, YlenRatios, YcellRatios, YexpRatios);
```

`splitEdge` mimics OpenFOAM multi-grading:

- `YL` is the total wall-normal length used to build the full `Y` vector.
- `YlenRatios` controls the relative physical lengths of the grading regions.
- `YcellRatios` controls how many cells each grading region receives.
- `YexpRatios` is the expansion ratio of the last cell to the first cell in
  each region.

Only `Y(1:NI+1)` is written into the OpenFOAM `points` file. Because
`NN = NI + NI_ADD`, the last `NI_ADD` wall-normal cells are available for
`suboff_mesh_2d_ADD.plt` but are not part of the OpenFOAM point field.

Take care with this distinction when changing `YL`, `NI`, or `NI_ADD`. If the
far-field location in OpenFOAM is not where expected, check the actual maximum
of `Y(1:NI+1)`, not only `YL`.

## Streamwise Surface Spacing

The streamwise spacing is built along hull arc length:

```matlab
SL = 0.65;
NS = NL/2 + NJ + NJ_ADD;

SlenRatios  = [0.35 , 0.95 , 9 , 12];
ScellRatios = [NL/2 , NJ1 , NJ2, NJ3 + NJ_ADD];
SexpRatios  = [1, 1.1 , 4 , 1];

S = splitEdge(SL, NS, SlenRatios, ScellRatios, SexpRatios);
X = Get_x_from_S(S);
```

`S` is distance along the meridional hull curve. `Get_x_from_S` converts arc
length back to axial coordinate:

- On the bow, it solves `S(x) = S_target` using `fzero`.
- On the parallel mid-body, where `dR/dx = 0`, it uses the exact relation
  `x = x_bow_end + (S - S_bow_end)`.

`SL` therefore controls how far downstream the generated bow/forebody grid
extends. If `SL` is changed, the matching OpenFOAM topology must still agree
with the resulting point counts.

## SUBOFF Geometry Mathematics

The code uses the SUBOFF bow and parallel-body formula from the appendix shown
in the supplied image. In the original nondimensional foot-based coordinates,

```text
R0(x) = Rmax * U(x)^alpha
Rmax = 5/6
alpha = 1/2.1
A = 1.126395101
B = 0.442874707
```

where, for the bow region `0 <= x <= 10/3`,

```text
U(x) = A*x*(0.3*x - 1)^4
     + B*x^2*(0.3*x - 1)^3
     + 1
     - (0.3*x - 1)^4*(1.2*x + 1)
```

For `x > 10/3`, the radius is constant:

```text
R0(x) = Rmax
```

The MATLAB implementation scales this geometry to meters:

```matlab
RATIO = 0.2;
x_dimless = x / RATIO / 0.3048;
R = R0 * 0.3048 * RATIO;
```

With the default `RATIO = 0.2`:

- Bow end: `(10/3) * 0.3048 * 0.2 = 0.2032 m`
- Maximum radius: `(5/6) * 0.3048 * 0.2 = 0.0508 m`

The slope used for normal extrusion is:

```text
K = dR/dx
```

and the arc length is:

```text
S(x) = integral sqrt(1 + K(x)^2) dx
```

## Point Extrusion Mathematics

For a surface point `(XH, YH, ZH)`, the local radius is:

```text
R = sqrt(YH^2 + ZH^2)
```

Given wall-normal distance `DY` and slope `K = dR/dx`, `Get_P.m` computes the
offset point:

```text
XP = XH - DY*K / sqrt(1 + K^2)
YP = (R + DY / sqrt(1 + K^2)) * YH/R
ZP = (R + DY / sqrt(1 + K^2)) * ZH/R
```

This is the outward normal in the meridional `(x, r)` plane, then revolved back
into three dimensions using the original circumferential direction. For the
special centerline point of the square H-block, the code writes:

```matlab
XP = -Y(I);
YP = 0.0;
ZP = 0.0;
```

This creates the upstream centerline branch needed by the H-type bow topology.

## Surface Construction Details

### Rectangular H-Block Bow Cap

`RECT_Y` and `RECT_Z` are square coordinates spanning:

```matlab
RECT_L = 2 * Get_R(X(NL/2+1));
```

The square side length is tied to the body radius at the H-O transition. The
corresponding `RECT_X` is obtained from `Get_x_from_R(RECT_R)`, which maps each
square radius back to the bow axial position. This works because the SUBOFF bow
radius is monotonic from the tip to the parallel body.

### Cylindrical/O-Type Surface

The cylindrical surface is assembled in three streamwise sections:

- `CY1`: transition from rectangular H-block perimeter to circular O-grid.
- `CY2`: circular O-grid over the next hull segment.
- `CY3`: circular O-grid over the final hull segment.

`CY1` performs the angular blending from the square corner angles to the
circular angles. `CY2` and `CY3` use a direct circumferential sweep:

```matlab
theta = linspace(-pi/4, 7*pi/4, NK+1);
Y = R(x) cos(theta)
Z = R(x) sin(theta)
```

The final arrays are concatenated as:

```matlab
CY_X = [CY_X1(:,1:NJ1) CY_X2(:,1:NJ2) CY_X3];
CY_Y = [CY_Y1(:,1:NJ1) CY_Y2(:,1:NJ2) CY_Y3];
CY_Z = [CY_Z1(:,1:NJ1) CY_Z2(:,1:NJ2) CY_Z3];
```

The duplicate columns are omitted at internal section boundaries to avoid
repeated streamwise planes.

## Output Files

Running `grid.m` writes the following files.

| File | Meaning |
| --- | --- |
| `points` | OpenFOAM v2312 `vectorField` containing all generated point coordinates. This is the file intended for `constant/polyMesh/points`. |
| `Tecplot_InputFiles/Y.txt` | First `NI+1` wall-normal distances written for checking the OpenFOAM grid spacing. |
| `Tecplot_InputFiles/RECT_surface.plt` | Tecplot surface file for the square H-type surface patch, including Tecplot headers. |
| `Tecplot_InputFiles/RECT_surface_XYZ.plt` | Raw XYZ version of the square H-type surface patch. |
| `Tecplot_InputFiles/CY1_surface.plt` | Tecplot surface file for the H-to-O transition section. |
| `Tecplot_InputFiles/CY2_surface.plt` | Tecplot surface file for the second O-type surface section. |
| `Tecplot_InputFiles/CY3_surface.plt` | Tecplot surface file for the third O-type surface section. |
| `Tecplot_InputFiles/CY_surface.plt` | Combined cylindrical/O-type surface with Tecplot headers. |
| `Tecplot_InputFiles/CY_surface_XYZ.plt` | Raw XYZ version of the combined cylindrical/O-type surface. |
| `Tecplot_InputFiles/suboff_mesh_2d.plt` | 2D meridional grid section using only the `NI+1` wall-normal nodes used by OpenFOAM. |
| `Tecplot_InputFiles/suboff_mesh_2d_ADD.plt` | 2D meridional grid section including the extra `NI_ADD` wall-normal nodes and `NJ_ADD` streamwise node. Used for 2d calculation in the extended grading. |

## OpenFOAM Point Ordering

The OpenFOAM `points` file is written in two blocks.

First, the O-type/cylindrical surface layers:

```matlab
for K = 1:(NL*4)       % same as NK circumferential stations
    for J = 1:(NJ+1)   % streamwise surface stations
        for I = NI+1:-1:1
            ...
        end
    end
end
```

Second, the interior square H-block points:

```matlab
for J = NL:-1:2
    for K = 2:NL
        for I = NI+1:-1:1
            ...
        end
    end
end
```

The square block writes only the interior `(NL-1)^2` points per wall-normal
layer because the H-block perimeter is already represented by the surrounding
O-type block.

If you change loop order, dimensions, or indexing, the OpenFOAM connectivity
files must be regenerated to match. A valid `points` file alone is not enough
to create a valid OpenFOAM mesh.

## Adapting to Another Axisymmetric Body

To use this generator for another body of revolution, keep the topology and
spacing logic but replace the geometry mappings:

- `Get_R(x)`: return body radius at axial coordinate `x`.
- `Get_K(x)`: return `dR/dx` in the same dimensional units as `x` and `R`.
- `Get_S(x)`: return meridional arc length from the nose to `x`.
- `Get_x_from_S(S)`: invert the arc-length mapping.
- `Get_x_from_R(R)`: invert the radius mapping over the bow region used by the
  rectangular H-block. This requires the bow radius to be one-to-one over the
  H-block mapping interval.

The current H-block construction assumes the bow radius increases monotonically
from zero to the transition radius. For a bulbous bow or a body with a
non-monotonic nose radius, `Get_x_from_R` and the rectangular cap construction
must be redesigned.

## Parameters That Need Extra Care

- Keep `NK` divisible by 8. The code uses `NL = NK/4` and `NL/2 = NK/8` in
  array sizes and sector mirroring.
- Keep `NL = NK/4` unless all dependent loops and point-count formulas are
  updated.
- Keep `NI`, `NJ1`, `NJ2`, `NJ3`, `NJ_ADD`, and `NI_ADD` consistent with any
  OpenFOAM connectivity files.
- Check `P_NUM` before running at high resolution. ASCII OpenFOAM point files
  become very large.
- Make sure `RATIO` is changed consistently in `Get_R.m`, `Get_K.m`,
  `Get_S.m`, `Get_x_from_S.m`, and `Get_x_from_R.m`.
- Do not use `Get_P.m` at a point where `R = 0`, except for the explicit
  centerline special case already handled in `grid.m`.
- Confirm the transition point `X(NL/2+1)` lies on the monotonic bow section.
  `RECT_X = Get_x_from_R(RECT_R)` depends on that.
- When changing grading ratios, inspect `suboff_mesh_2d.plt` and
  `suboff_mesh_2d_ADD.plt` before using the OpenFOAM mesh.

## Suggested Validation Workflow

Before using the output in a CFD case:

1. Run a low-resolution test mesh.
2. Open `suboff_mesh_2d.plt` in Tecplot or ParaView-compatible tooling and
   check wall-normal spacing and streamwise clustering.
3. Open `RECT_surface.plt` and verify that the square H-block cap is smooth and
   centered on the nose.
4. Open `CY_surface.plt` and check that the circumference closes cleanly.
5. Check the first and last circumferential planes for duplicate or missing
   points.
6. Copy `points` into `constant/polyMesh/points` only after confirming that the
   rest of the OpenFOAM `polyMesh` files were generated for the same dimensions
   and ordering.

## Known Limitations

- The project currently writes only the OpenFOAM `points` file.
- Geometry constants are duplicated across helper functions, so changing the
  body scale or shape requires coordinated edits.
- The generated OpenFOAM file is ASCII, which can be extremely large for the
  default resolution.
- The code is specialized to a monotonic SUBOFF-style bow followed by a
  parallel mid-body. Other bodies may require changes to the inverse geometry
  functions.
