# 08_massflow

This folder corrects the turbulent inlet velocity time series so the inlet mass
flow is consistent at every time step. It uses the mass-flow history measured
from `07_mesh_for_turb`, blends the synthetic turbulence from `06_tg` with the
non-uniform mean inlet from `05_ini_2d_RAS`, rescales each time slice to a
reference flux, and writes OpenFOAM `constant/boundaryData/INLET/<time>/U`
files.

The main downstream product is a corrected `constant/boundaryData/INLET`
directory containing `points` and time-resolved `U` files for
`timeVaryingMappedFixedValue`.

## Workflow Idea

`06_tg` generates a turbulent velocity sequence, but the instantaneous inlet
flux can drift. `07_mesh_for_turb` measures that flux with
`system/mass_flow_cal` and writes `postProcessing/inletletMassFlow/...`.
This folder uses that measured `sum(phi)` history to correct the inlet sequence.

For each inlet point, `inlet_code/suboff_main.F90` applies a radial weight
function:

```text
U_mix = (U_tg - Uinf e_x) * weight + U_RECT
```

where `U_tg` is the synthetic turbulent inlet from `06_tg/vel`, `U_RECT` is the
non-uniform mean field from `PU_RECT_SUR.plt`, and `weight` tapers the
turbulent fluctuation near the outer radius. The field is then scaled by:

```text
U_corrected = U_mix * PHI_REF / abs(phi(t))
```

The final `MID` time steps are linearly smoothed back toward the first corrected
time step so the sequence can be repeated more smoothly.

## Repository Layout

```text
08_massflow/
`-- constant/
    `-- boundaryData/
        |-- INLET/                    # Corrected OpenFOAM mapped boundary data
        |-- phi_ref/
        |   |-- ReadData/
        |   |   `-- surfaceFieldValue.dat
        |   |-- suboff_main.F90       # Computes mean reference mass flow
        |   `-- run.sh
        |-- mkdir_code/
        |   |-- mkdir.F90             # Creates INLET time directories
        |   |-- parameter.h
        |   `-- run.sh
        |-- points_code/
        |   |-- cylinder.F90          # Writes INLET/points from INLET1.plt
        |   |-- parameter.h
        |   |-- run.sh
        |   `-- ReadData/
        `-- inlet_code/
            |-- suboff_main.F90       # Corrects and writes INLET/<time>/U
            |-- parameter.h
            |-- param.fi
            |-- head.fi
            |-- common.fi
            |-- run.sh
            |-- ofrun.sh
            |-- ReadData/
            `-- Tecplot_InputFiles/
```

## Quick Start

Run each command from the directory named in the step.

1. Prepare the mass-flow reference.

   Copy the measured mass-flow file from `07_mesh_for_turb`:

   ```text
   07_mesh_for_turb/postProcessing/inletletMassFlow/0/surfaceFieldValue.dat
   ```

   into:

   ```text
   constant/boundaryData/phi_ref/ReadData/surfaceFieldValue.dat
   constant/boundaryData/inlet_code/ReadData/surfaceFieldValue.dat
   ```

   Then compute the time-mean flux:

   ```bash
   cd constant/boundaryData/phi_ref
   ./run.sh
   ```

   This writes:

   ```text
   phi_ref.txt
   ```

   Note: `inlet_code/suboff_main.F90` currently hard-codes
   `PHI_REF=5.623942839039030E-002`, so update the source if you want to use the
   newly computed `phi_ref.txt` value.

2. Create the `INLET` time directories.

   ```bash
   cd constant/boundaryData/mkdir_code
   ./run.sh
   ```

   `mkdir.F90` generates an `ofrun.sh` with `mkdir -p ../INLET/<time>` commands
   for `0` through `DT*NT`, then executes it.

3. Write the mapped inlet points.

   Copy:

   ```text
   03_get_inlet1_points/INLET1.plt
   ```

   into:

   ```text
   constant/boundaryData/points_code/ReadData/INLET1.plt
   ```

   Then run:

   ```bash
   cd constant/boundaryData/points_code
   ./run.sh
   ```

   This writes:

   ```text
   constant/boundaryData/INLET/points
   ```

4. Prepare the corrected-inlet inputs.

   Copy:

   ```text
   05_ini_2d_RAS/get_ini_p_u/PU_RECT_SUR.plt
   ```

   into:

   ```text
   constant/boundaryData/inlet_code/ReadData/PU_RECT_SUR.plt
   ```

   Confirm that the turbulent velocity files exist at the relative path used by
   `inlet_code/suboff_main.F90`:

   ```text
   ../../../../06_tg/vel/u000001.plt
   ../../../../06_tg/vel/u000002.plt
   ...
   ```

5. Generate the corrected OpenFOAM inlet data.

   ```bash
   cd constant/boundaryData/inlet_code
   ./run.sh
   ```

   `run.sh` builds `a.zf` and submits `ofrun.sh`, which runs:

   ```bash
   mpirun -np 100 ./a.zf
   ```

   The output files are written under:

   ```text
   constant/boundaryData/INLET/<time>/U
   ```

6. Optional: create Tecplot preplot scripts for inlet animation.

   ```bash
   cd constant/boundaryData/inlet_code/Tecplot_InputFiles/UINLET_Preplot
   ./run.sh
   ```

   This regenerates `run_preplot.ps1`, which converts
   `../UINLET/UINLET_Animate_*.plt` files with Tecplot `preplot`.

## Input Files

| File | Used By | Meaning |
| --- | --- | --- |
| `phi_ref/ReadData/surfaceFieldValue.dat` | `phi_ref/suboff_main.F90` | Mass-flow history copied from `07_mesh_for_turb/postProcessing/inletletMassFlow/0/`. |
| `inlet_code/ReadData/surfaceFieldValue.dat` | `inlet_code/suboff_main.F90` | Per-time-step `sum(phi)` used to rescale the inlet velocity field. |
| `inlet_code/ReadData/PU_RECT_SUR.plt` | `inlet_code/suboff_main.F90` | Non-uniform mean pressure/velocity field on `INLET1`, generated by `05_ini_2d_RAS`. |
| `points_code/ReadData/INLET1.plt` | `points_code/cylinder.F90` | Inlet point coordinates from `03_get_inlet1_points`. |
| `06_tg/vel/u######.plt` | `inlet_code/suboff_main.F90` | Synthetic turbulent inlet velocity sequence. |

## Output Files

| File | Meaning |
| --- | --- |
| `phi_ref/phi_ref.txt` | Mean `sum(phi)` computed from `surfaceFieldValue.dat`; informational unless copied into the hard-coded `PHI_REF`. |
| `constant/boundaryData/INLET/points` | OpenFOAM mapped-boundary point coordinates. |
| `constant/boundaryData/INLET/0/U` | Initial inlet field, written from `PU_RECT_SUR.plt` without turbulent fluctuation. |
| `constant/boundaryData/INLET/<time>/U` | Mass-flow-corrected turbulent inlet velocity field for each time directory. |
| `inlet_code/RATIO.plt` | Raw Tecplot-style radial weight function output. |
| `inlet_code/Tecplot_InputFiles/ITIweightFunc.plt` | Tecplot visualization of the radial blending weight. |
| `inlet_code/Tecplot_InputFiles/UINLET/UINLET_Animate_<rank>.plt` | Tecplot animation slices written by each MPI rank for inspection. |
| `inlet_code/Tecplot_InputFiles/UINLET_Preplot/run_preplot.ps1` | PowerShell script for Tecplot `preplot` conversion. |

## Important Parameters

| Location | Parameter | Meaning |
| --- | --- | --- |
| `mkdir_code/parameter.h` | `DT=2e-4`, `NT=5000` | Time directory spacing and number of steps. |
| `inlet_code/parameter.h` | `DT=2e-4`, `NT=5000`, `ST=0`, `MID=20` | Velocity-file mapping and final smoothing window. |
| `inlet_code/parameter.h` | `RMAX=0.09`, `DELTA_RATIO=0.3` | Radius and transition width for the fluctuation blending weight. |
| `inlet_code/parameter.h` | `N2=160`, `N3=160`, `NALL=N2*N3` | Square inlet dimensions; must match `INLET1` and `06_tg`. |
| `inlet_code/param.fi` | `NPR=100`, `NPC=1`, `NPP=100` | MPI process count for corrected-inlet generation. |
| `points_code/parameter.h` | `N=160`, `NALL=N*N` | Mapped boundary point count. |
| `inlet_code/suboff_main.F90` | `UINF=1.649194` | Uniform velocity subtracted from `06_tg` before blending. |
| `inlet_code/suboff_main.F90` | `PHI_REF=5.623942839039030E-002` | Target reference flux used to rescale every time step. |

Keep `DT` synchronized with `06_tg/parameter.h`, `07_mesh_for_turb/controlDict`,
and any downstream OpenFOAM case that reads the corrected `boundaryData`.

## Correction Details

The radial weight function is based on the inlet-point radius `R`:

```text
TEMP1 = DELTA_RATIO * 0.5 * RMAX
TEMP2 = R - (RMAX - TEMP1)
```

The weight is:

```text
1,                                      TEMP2 <= -TEMP1
0,                                      TEMP2 >=  TEMP1
1 - 0.5*(1 + TEMP2/TEMP1 + sin(pi*TEMP2/TEMP1)/pi), otherwise
```

For most time steps, the corrected inlet is:

```text
Ux = ((Utg_x - UINF) * weight + Urect_x) * PHI_REF / abs(phi(t))
Uy = ( Utg_y         * weight + Urect_y) * PHI_REF / abs(phi(t))
Uz = ( Utg_z         * weight + Urect_z) * PHI_REF / abs(phi(t))
```

For the last `MID` time steps, the code linearly blends from the end field back
to the first corrected field before applying the flux scaling. This makes the
generated sequence smoother if it is repeated periodically.

## BoundaryData Notes

The generated directory is intended to be consumed by an OpenFOAM boundary
condition like:

```text
type timeVaryingMappedFixedValue;
```

The folder name must match the patch name exactly:

```text
constant/boundaryData/INLET
```

and the `points` file and all time-directory `U` files must use the same point
ordering as `INLET1.plt`.

## Notes And Limitations

- No `Makefile` files were found in this folder snapshot, but the `run.sh`
  scripts call `make`. Restore or provide the local build rules before running
  the helper programs.
- `inlet_code/suboff_main.F90` uses fixed relative paths to `06_tg/vel` and
  fixed file names under `ReadData`; keep the workflow layout unchanged or
  update those paths.
- `PHI_REF` is hard-coded in `inlet_code/suboff_main.F90`; recomputing
  `phi_ref.txt` does not automatically change the correction target.
- `NT` should match the number of available `06_tg/vel/u######.plt` files.
- The MPI converter assumes `NT/NPP` work per rank. Keep `NT`, `NPP`, and the
  `mpirun -np` value in `ofrun.sh` consistent.
