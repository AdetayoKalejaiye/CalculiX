# CalculiX Relational Map

This map describes the principal relationships in the repository. It is
intentionally organized around runtime data flow and build-time dependencies,
because the implementation is a large flat C/Fortran source tree rather than
a package hierarchy.

## Runtime and data-flow map

```mermaid
flowchart TD
    CLI["CalculiX executable<br/>src/CalculiX.c: main"] --> OPEN["Open job files"]
    CLI --> READ["readinput.c<br/>parse .inp and build pointers"]
    READ --> INIT["ini_cal.c<br/>initialize counters and flags"]
    INIT --> ALLOC["allocation.f<br/>allocate shared model/state arrays"]
    ALLOC --> STEP["calinput.f<br/>interpret analysis steps"]

    STEP --> MODEL["Mesh, sets, materials,<br/>loads, BCs, MPCs"]
    STEP --> CONTACT["Contact and tied/mortar<br/>constraints"]
    STEP --> COUPLED["Thermal, fluid, CFD,<br/>network, electromagnetic data"]

    MODEL --> ASSEMBLY["Element and matrix assembly<br/>mafills*, matrixassembles.f"]
    CONTACT --> ASSEMBLY
    COUPLED --> ASSEMBLY

    ASSEMBLY --> STRUCT["Sparse structure/storage<br/>mastruct*, matrixstorage*"]
    STRUCT --> SOLVE{"Analysis / solver path"}

    SOLVE --> STATIC["Static/nonlinear<br/>statics.f"]
    SOLVE --> DYNAMIC["Dynamic/transient<br/>dynamics.f, dyna.c"]
    SOLVE --> MODAL["Modal/eigenvalue<br/>modaldynamics.f, arpack*.c"]
    SOLVE --> MULTIPHYS["Thermal/fluid/network/<br/>electromagnetic drivers"]

    STATIC --> BACKEND["Equation backend"]
    DYNAMIC --> BACKEND
    MODAL --> BACKEND
    MULTIPHYS --> BACKEND

    BACKEND --> ITER["solveeq.c / iterative routines"]
    BACKEND --> SPOOLES["spooles.c"]
    BACKEND --> PARDISO["pardiso.c"]
    BACKEND --> PASTIX["pastix.c"]

    BACKEND --> UPDATE["Convergence, time-step,<br/>and state update"]
    UPDATE --> STEP
    UPDATE --> RESULTS["results*.c/f<br/>recover physical results"]
    RESULTS --> DAT[".dat writers"]
    RESULTS --> FRD["frd*.c/f<br/>.frd visualization output"]
    UPDATE --> RESTART["restart read/write"]
```

## Build and dependency map

```mermaid
flowchart LR
    CMAKE["CMakeLists.txt"] --> SOURCES["src/*.c + src/*.f"]
    PRESETS["CMakePresets.json"] --> CMAKE
    FINDERS["cmake/Find*.cmake"] --> CMAKE
    COMPAT["cmake/compat"] --> CMAKE

    SOURCES --> CCX["ccx library"]
    ENTRY["src/CalculiX.c"] --> EXE["CalculiX executable"]
    CCX --> EXE

    OPENMP["OpenMP"] --> CCX
    BLAS["BLAS/LAPACK or Intel MKL"] --> CCX
    ARPACK["ARPACK / ARPACK-ng"] --> CCX
    SPOOLES["SPOOLES"] --> CCX
    PARDISO["PARDISO via MKL"] --> CCX
    PASTIX["PaStiX or PaStiX4CalculiX"] --> CCX

    ILP64["CALCULIX_INT64 / INTSIZE64"] -. changes ITG ABI .-> SOURCES
    ILP64 -. requires matching integer interfaces .-> BLAS
    ILP64 -. requires matching integer interfaces .-> ARPACK
    ILP64 -. requires matching integer interfaces .-> PASTIX

    LEGACY["src/Makefile + Makefile.inc"] --> SOURCES
    LEGACY --> EXTERNAL["Manual sibling dependencies"]
```

## Repository relationship map

```mermaid
flowchart TD
    ROOT["Repository"] --> CORE["src/<br/>solver implementation"]
    ROOT --> BUILD["CMakeLists.txt<br/>CMakePresets.json<br/>cmake/"]
    ROOT --> TESTS["test/<br/>input decks and references"]
    ROOT --> MANUAL["doc/<br/>LaTeX manual and figures"]
    ROOT --> INSTALL["INSTALL.md<br/>build/dependency instructions"]

    BUILD --> CORE
    BUILD --> DEPS["Compiler, OpenMP,<br/>BLAS/LAPACK, ARPACK,<br/>sparse solvers"]
    CORE --> OUTPUT[".dat, .frd, restart,<br/>matrix and feature outputs"]
    TESTS --> CORE
    CORE --> TESTS
    MANUAL -. explains user-visible behavior .-> CORE
    INSTALL -. explains how to build .-> BUILD
```

## Key relationship notes

| Relationship | Meaning |
| --- | --- |
| `CalculiX.c` → shared arrays | The executable owns most model and solution state and passes it through broad C/Fortran interfaces. |
| `CalculiX.h` ↔ Fortran routines | The header defines naming, integer-size, allocation, and function-signature conventions for the language boundary. |
| `calinput.f` → analysis drivers | Input keywords select and configure static, dynamic, modal, coupled-physics, and specialized procedures. |
| Assembly → solver adapter | Matrix structure and values are prepared once per required phase, then passed to an iterative or external sparse backend. |
| Solver/update → results | Converged state is postprocessed into physical quantities before file writers emit user-facing results. |
| `CMakeLists.txt` → optional adapters | Build flags determine which external solver implementation is compiled and linked. |
| `test/*.inp` ↔ `.ref` files | Regression cases define expected numerical and output behavior for the executable. |
| `CalculiXstep.c` ↔ core routines | The embedding-style API shares the same numerical machinery but preserves state across caller-driven step invocations. |
