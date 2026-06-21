# Project Experience Notes

## Cursor Working Rules

- Read this file at the beginning of every new chat before planning or editing.
- Only modify files inside ``d:\OpenFOAM\v2312\SUBOFF` unless the user explicitly
  approves another exact path. 
- Preserve user edits and generated case state. Do not reset, delete processor
  folders, clean meshes, or overwrite OpenFOAM time directories unless the user
  has asked for that operation.
- When editing workflow logic, update the nearby README or this experience file
  with any new assumptions about paths, dimensions, solver settings, or output
  handoff files.

## Project Workflow Notes

This directory combines several kinds of files:

- MATLAB scripts generate structured grid coordinates and post-process boundary
  layer, curvature, FIK, and mesh data.
- Fortran/MPI programs generate or transform boundary data, point data, and
  Tecplot/OpenFOAM input files.
- OpenFOAM v2312 case directories contain `system`, `constant`, `0`/`0.orig`,
  `constant/polyMesh`, and `constant/boundaryData` workflows.
- Shell scripts compile tools, run OpenFOAM utilities, decompose/reconstruct
  cases, submit jobs, and clean selected generated directories.

Many files depend on matching dimensions and ordering. Before changing
resolution, array loops, point ordering, boundary names, or processor counts,
check the related MATLAB, Fortran, OpenFOAM dictionary, and README notes
together.

## MATLAB Coding Notes

- MATLAB mesh scripts usually rely on the current working directory. If a script
  opens `Tecplot_InputFiles/...`, `points`, `Y.txt`, or similar relative paths,
  document and preserve the expected run directory.
- Keep derived mesh dimensions close to their source variables. Common examples
  are `NL = NK/4`, `NJ = NJ1 + NJ2 + NJ3`, `NN = NI + NI_ADD`, and
  `NS = NL/2 + NJ + NJ_ADD`.
- Preserve divisibility and indexing assumptions. In the H-O grid generator,
  `NK` must remain divisible by 8 because the code uses `NK/4` and `NK/8`
  sectors for square-to-circular mapping.
- Distinguish OpenFOAM nodes from extra visualization or inflow-helper nodes.
  For example, `NI_ADD` and `NJ_ADD` may extend 2D/Tecplot helper output without
  being written into the final OpenFOAM `points` field.
- Geometry constants such as `RATIO`, SUBOFF bow coefficients, radius formulas,
  slope `Get_K`, arc length `Get_S`, and inverse mappings must be updated
  consistently across helper functions.
- Use structured numeric output formats already present in the scripts, such as
  `%22.15E`, unless a downstream reader is changed too.
- For Tecplot files, keep headers, `VARIABLES`, `ZONE`, and duplicated closing
  circumferential points consistent with the existing plotting workflow.

## Fortran Coding Notes

- Most Fortran tools are fixed around included configuration files such as
  `head.fi` and `parameter.h`. Check these include files before changing array
  sizes, MPI process counts, or output dimensions.
- Preserve MPI topology assumptions. Programs may use `MPI_CART_CREATE`,
  `MPI_CART_COORDS`, `MPI_CART_SUB`, and process-count checks such as
  `NPP == iMPI_NumProcs`; these must agree with run scripts and job settings.
- Many readers/writers assume exact nested loop order, often `K`, `J`, then
  component or point values. Changing loop order can break Tecplot files,
  OpenFOAM boundary data, or MATLAB post-processing.
- Be careful with large fixed-length character buffers and formatted writes.
  Some files emit OpenFOAM list syntax or long boundary records, so trimming or
  shortening buffers can corrupt output.
- Keep output responsibility clear: raw files such as `PU_RECT_SUR.plt`,
  `RATIO.plt`, or OpenFOAM `boundaryData` lists may be consumed later, while
  files under `Tecplot_InputFiles` are usually inspection outputs.
- Existing code contains English and Chinese comments. Preserve meaningful
  comments; add concise comments only when they clarify non-obvious indexing,
  MPI decomposition, or file format assumptions.

## OpenFOAM v2312 Notes

- Keep OpenFOAM dictionary headers consistent with v2312 style and the correct
  `FoamFile` class/object names. Do not change dictionary object names unless
  the file is being renamed or the case structure changes.
- Existing cases commonly use ASCII output with high precision, for example
  `writeFormat ascii` and `writePrecision 14`. Generated mesh and boundary-data
  files should remain compatible with OpenFOAM v2312 readers.
- `decomposeParDict` processor counts must match run scripts. For example, a
  case using `numberOfSubdomains 384` must agree with `mpirun -np 384`,
  processor-folder cleanup loops, and any Fortran/MPI decomposition constants.
- Mesh-construction utilities such as `redistributePar`, `reconstructPar`, 
  `foamCleanPolyMesh`, and manual `processor*` folder moves can overwrite or remove case state. Treat them as operational commands, not harmless formatting checks.
- `blockMeshDict`, and `constant/polyMesh/*` files depend on the
  same geometry and indexing as generated `points`, `faces`, `owner`,
  `neighbour`, and `boundary`. A valid `points` file alone is not enough.
- Boundary names and folder names in `constant/boundaryData` must match
  OpenFOAM patch names exactly. Preserve capitalization and directory structure.
- When changing solver dictionaries, note whether the case is laminar, RAS/RANS,
  transient, or steady. Preserve the intended solver such as `simpleFoam` or
  `pimpleFoam` unless the workflow is being intentionally changed.

## Shell And Run Script Notes

- Run scripts assume Bash syntax even on Windows-hosted workspaces, including
  brace loops like `{0..15}` and OpenFOAM environment commands. Keep scripts
  portable to the actual OpenFOAM shell environment, not PowerShell.
- Before editing processor loops, update every related `-np`, processor-folder
  loop, decomposition dictionary, and MPI/Fortran constant.
- Treat `rm -rf`, mesh cleaning, reconstruction, and case-copy commands as
  destructive. Document what state they remove or recreate.
- If a script compiles a Fortran tool and then checks for an output executable,
  keep the failure check before job submission or OpenFOAM execution.
- Prefer explicit stdout log files such as `log.run` for long solver runs, and
  avoid hiding failures by chaining too many critical commands without checks.

## README Writing Notes

A README should explain not only what a folder contains, but also
how its files fit into the SUBOFF/OpenFOAM mesh-construction workflow.

### Recommended Structure

Start with a clear title naming the folder or tool, followed by a short
paragraph that states its purpose and main output. Mention the target format
early, for example OpenFOAM `constant/polyMesh/*` files, Tecplot `.plt` files,
MATLAB scripts, or OpenFOAM case folders.

For technical folders, use sections similar to:

- `Grid Idea`, `Case Idea`, or `Workflow Idea`: explain the method/topology in
  plain language before listing files.
- `Repository Layout`: show the folder tree and briefly describe each important
  script, case directory, or generated output.
- `Quick Start`: provide exact steps and the expected working directory. Relative
  paths matter in this project, so state where commands or scripts must be run.
- `Important Variables` or `Important Parameters`: list key resolution,
  geometry, grading, solver, or parallel-decomposition controls.
- `Output Files`: separate visualization/check files from downstream files used
  by later OpenFOAM steps.
- `Ordering`, `Indexing`, or `Connectivity`: document any point, face, cell, or
  processor ordering that another file depends on.
- `Parameters That Need Extra Care`: collect constraints that can silently break
  the mesh or solver.
- `Known Limitations`: state what the folder does not yet support.

### Content Rules

- Always document the relation between generated files and the next workflow
  step. For example, say whether a file is only for Tecplot checking or is meant
  to be copied into `constant/polyMesh`.
- When a script writes files using relative paths, explicitly state the required
  run directory and where the outputs will appear.
- Include formulas for point counts, cell counts, or derived dimensions when
  they determine compatibility with OpenFOAM topology files.
- Explain derived variables near their source variables. For example, if
  `NL = NK/4`, document both the meaning and the required divisibility of `NK`.
- Record invariants, not just default values. Good README notes say which values
  may be changed freely and which require regenerating related files.
- Warn when multiple files duplicate the same geometry or scale constants, since
  inconsistent edits are easy to make.
- Keep MATLAB/OpenFOAM command snippets copyable and minimal. Avoid long code
  dumps unless they show a loop order, formula, or command sequence that users
  must preserve.
- Prefer tables for variable meanings and output-file meanings when the list is
  long enough to compare entries.
- Use consistent file names and capitalization exactly as they appear on disk,
  including names such as `suboff_mesh_2d_ADD.plt`.

### Style Notes

- Write for someone rerunning or modifying the case later, not only for someone
  reading the code today.
- Put the practical warnings close to the operation they affect, then repeat the
  most important constraints in a final care/limitations section.
- Use concise paragraphs before bullet lists so readers understand the purpose
  before scanning details.
- Prefer concrete project terms such as H-block, O-grid, wall-normal spacing,
  streamwise spacing, `points`, `faces`, `owner`, `neighbour`, and `boundary`.
- If a section is unfinished, mark it clearly instead of leaving vague text such
  as "Waiting for completed ...".

### Checklist Before Finishing A README

- The folder purpose and main output are clear in the first few lines.
- The user knows exactly where to run the script or OpenFOAM commands.
- Generated files are classified as visualization/check outputs or downstream
  workflow inputs.
- Key parameters, derived values, and divisibility/indexing constraints are
  documented.
- Any dependency on matching OpenFOAM connectivity or case files is stated.
- Known limitations and risky edits are called out explicitly.
