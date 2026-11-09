# CalculiX Codebase Guide

This guide is a map for developers working on the CalculiX CrunchiX solver in
this repository. CalculiX is a large mixed C/Fortran finite-element program:
C provides the executable boundary, shared conventions, memory helpers,
parallel/runtime glue, file-format code, and several solver adapters; Fortran
contains most input interpretation, element formulations, constitutive models,
analysis procedures, and numerical kernels.

For the visual overview, see [`relational-map.md`](relational-map.md). For
source-grounded explanations of the numerical procedures, see
[`algorithm-guide.md`](algorithm-guide.md).

## 1. Repository at a glance

| Area | Purpose |
| --- | --- |
| `src/` | Core solver. Nearly all production code lives in this single source directory. |
| `cmake/` | CMake flags, dependency finders, and compatibility adapters. |
| `CMakeLists.txt` | Modern build graph: creates the `ccx` library and `CalculiX` executable. |
| `CMakePresets.json` | GNU/Intel LLVM presets for LP64 and ILP64 builds. |
| `src/Makefile*` | Legacy/manual build variants, including multithreaded, integer-64, and MFront builds. |
| `test/` | Input decks and reference `.dat`/`.frd` outputs, plus comparison scripts. |
| `doc/` | LaTeX manual and figures. |
| `pictures/` | Images used by the top-level README. |
| `INSTALL.md` | Dependency and build instructions. |
| `LOGBOOK` | Project history and development notes. |

The source directory is intentionally flat. Treat related filename families as
subsystems rather than expecting one directory per feature.

## 2. Build and dependency model

### Recommended CMake build

The standard GNU LP64 build uses 32-bit solver indices and SPOOLES:

```sh
cmake --preset gnu-lp64
cmake --build --preset gnu-lp64 -j
```

The executable is then `build/gnu-lp64/CalculiX`. The presets also provide
GNU/Intel LLVM and LP64/ILP64 combinations. ILP64 is intended for very large
models and requires compatible 64-bit BLAS, ARPACK, and sparse-solver
dependencies.

The CMake target relationship is:

```text
src/*.c + src/*.f
        │
        ▼
      ccx  ── external libraries (OpenMP, BLAS/LAPACK, ARPACK, sparse solver)
        │
        ▼
   CalculiX  (src/CalculiX.c)
```

`CMakeLists.txt` discovers C and Fortran sources recursively, excludes the C
entry point and a small set of support sources from the library, then links the
`CalculiX` executable against `ccx`. Optional compile-time integrations are
selected with `SPOOLES`, `PARDISO`, `PASTIX`, and `PASTIX4CALCULIX`.

### Legacy Make builds

`src/Makefile` requires manually configured sibling dependency trees and builds
the source lists from `src/Makefile.inc`. Other Makefiles cover specialized
variants. Prefer CMake for reproducible new builds; use Make when matching an
existing downstream or historical build is required.

### Important ABI conventions

`src/CalculiX.h` is the C/Fortran boundary:

- `FORTRAN(name, args)` expands compiler/platform-specific symbol names.
- `ITG` is `int` normally and `long long` when `INTSIZE64` is enabled.
- Allocation macros route through the project memory helpers.
- C declarations must match Fortran argument order and array layout exactly.

When changing a cross-language function, update both its implementation and
the declaration/call sites. A wrong scalar type or argument order can compile
successfully and still corrupt solver state.

## 3. Runtime lifecycle

The normal executable path is:

1. `src/CalculiX.c:main` parses the job name and opens the input deck.
2. `readinput.c` reads and normalizes the `.inp` deck, including input
   pointers, included files, sets, and keyword ordering.
3. `ini_cal.c` initializes counters, analysis flags, and state defaults.
4. `allocation.f` allocates arrays for the model, loads, constraints,
   materials, contact, thermal/network data, and output state.
5. `calinput.f` interprets the input deck one analysis step at a time.
6. Mesh, sets, materials, loads, boundary conditions, MPCs, contact, and
   coupled-physics data are assembled into shared arrays.
7. An analysis driver such as `statics.f`, `dynamics.f`, modal routines, or
   network/thermal/electromagnetic routines assembles equations.
8. Element and matrix routines form the system; `solveeq.c` or a selected
   sparse solver adapter solves it.
9. State variables are updated and convergence/time-step checks run.
10. `results*`, `frd*`, and `write*` routines emit `.dat`, `.frd`, restart,
    matrix, eigenvalue, and feature-specific output.

The loop is repeated for increments, iterations, and analysis steps. Most
working data is passed through large shared arrays owned by the top-level
driver, which explains the broad function signatures in `CalculiX.h`.

`CalculiXstep.c` exposes a second, embedding-oriented step driver. It keeps
solver state in static storage between calls and accepts mesh, loads, boundary
conditions, and state arrays from the caller. Changes to the main execution
path should be checked against this API-style path when the affected routines
are shared.

## 4. Source subsystem map

### Input, model setup, and state

- `readinput.c`, `readnewmesh.c`, `readfrd.c`: input and mesh/result readers.
- `calinput.f`, `calinput_rfn.f`: keyword and step interpretation.
- `ini_cal.c`, `allocation.f`, `allocation_rfn.f`: initialization and sizing.
- `create*.f`, `initialmeshs.f`, `modelchanges.f`: derived topology and model
  structures.
- `check*.c/f`, `inputerror.f`, `inputwarning.f`: validation and diagnostics.

### Element and matrix computation

- `mafills*.c`, `mafills*.f`, `matrixassembles.f`: element contributions and
  global assembly.
- `e_*.f`, `shape*.f`, `beam*.f`, `shell*.f`, and section routines: element
  formulations, interpolation, integration, and section behavior.
- `mastruct*.c`, `matrixsort.c`, `matrixstorage.c`, `matrixstorage.h`: sparse
  structure and storage.

### Analysis drivers

- `statics.f`: static and nonlinear equilibrium procedures.
- `dynamics.f`, `dyna.c`: transient and dynamic procedures.
- `modaldynamics.f`, `steadystate.c`, `arpack*.c`: modal, frequency, and
  eigenvalue paths.
- `heattransfers.f`, `therm*.f`, `cfds.f`, `network*.f`: thermal, fluid, CFD,
  and network coupling.
- `electromagnetics.c`, `electromagneticss.f`, and induction result routines:
  electromagnetic analyses.

### Materials and physics

- `materials.f` and `material*.f`: material data and constitutive dispatch.
- `plastic*.f`, `creep*.f`, `hyperelastics.f`, `damage*.f`, `umat*.f`: material
  behavior families.
- User/external material support is wired through specialized source files and
  the MFront Makefile variant.

### Constraints and contact

- `constraints.f`, `mpcs.f`, `createinterfacempcs.f`: constraint equations.
- `contact*.c/f`, `contactpairs.f`, `contactmortar.c`, `stressmortar.c`,
  `tiedcontact.c`: contact detection, mortar/tied contact, and contact forces.
- Boundary conditions and loads are represented in shared arrays and applied
  during assembly and solution updates.

### Solvers and external libraries

- `solveeq.c`: iterative/lumped-system equation solve path.
- `spooles.c`, `pardiso.c`, `pastix.c`: sparse direct solver adapters.
- `arpack.c`, `arpackbu.c`, `arpackcs.c`: eigenvalue adapters.
- `cmake/compat/`: source substitutions needed for supported external APIs.

### Results and file formats

- `results.c` and `results*.c/f`: recover strains, stresses, forces, energies,
  temperatures, and feature-specific results.
- `frd*.c/f`: `.frd` visualization output.
- `write*.c/f`: `.dat`, restart, mesh, matrix, vector, sensitivity, and
  feature-specific writers.

## 5. Where to make a change

| Change | Start here | Also inspect |
| --- | --- | --- |
| New input keyword or keyword behavior | `readinput.c`, `calinput.f` | `inputinfo.f`, validation, `test/*.inp` |
| New analysis procedure | `calinput.f`, relevant analysis driver | assembly, solver, results, restart behavior |
| Element formulation | matching `e_*.f`, shape/integration routines | `mafills*`, element tests, output recovery |
| Material law | `materials.f`, matching `material*.f`/`umat*.f` | state allocation, tangent/stress recovery, tests |
| Contact or constraint behavior | `contact*.c/f`, `constraints.f`, `mpcs.f` | assembly, convergence, contact tests |
| Solver integration | adapter source and `CMakeLists.txt` | integer ABI, linker options, solver-specific tests |
| Output field or file format | `results*`, `frd*`, `write*` | output labels and reference files |
| Cross-language interface | `CalculiX.h` and both callers/callees | `ITG`, `FORTRAN`, array dimensions, LP64/ILP64 |
| Build/dependency behavior | `CMakeLists.txt`, `CMakePresets.json`, `cmake/` | `INSTALL.md`, legacy Makefiles |

Avoid adding an isolated helper when a matching filename family already exists.
The existing code organizes behavior primarily by operation and physics
family, not by class or package.

## 6. Tests and validation

The `test/` directory is a regression corpus of input decks and expected
outputs. `test/README` documents the two main runners:

- `compare` runs serially.
- `compare_par` runs cases in parallel and compares results serially.

The comparison scripts check generated `.dat` and `.frd` files against
`.ref` files and create `error*` files for differences. Choose focused cases
that exercise the changed subsystem, then run the broader comparison suite when
the change affects shared assembly, solver, input, or output code.

Useful test families include beams and solids, contact/mortar, thermal and
fluid/network analyses, dynamics, electromagnetics, fracture, sensitivity,
and user-material cases. Keep reference outputs synchronized with intentional
behavior changes; do not update them merely to hide a regression.

## 7. Practical debugging strategy

1. Reproduce with the smallest existing `.inp` case that reaches the changed
   path.
2. Identify the analysis driver and the shared arrays it consumes and writes.
3. Check both LP64 and ILP64 assumptions if the change touches indices,
   sparse storage, or external libraries.
4. For numerical differences, distinguish input parsing, assembly, solver
   convergence, state update, and output recovery before changing tolerances.
5. For C/Fortran changes, verify symbol naming, argument order, dimensions,
   and allocation lifetime.
6. Compare `.dat` and `.frd` separately; a correct solver result can still
   expose an output-only regression.

## 8. Documentation and ownership boundaries

Use `INSTALL.md` for dependency/build details, `doc/CalculiX.tex` for the
user-facing theory and keyword manual, and this guide for repository structure
and code navigation. Keep architecture statements here tied to stable source
boundaries; detailed algorithm explanations belong next to the relevant
routine or in the manual.
