# CalculiX Algorithm Guide

This document explains the main algorithms implemented in this repository.
The descriptions are derived from the source code, comments, data structures,
and call sites; they are not a replacement for the CalculiX theory manual.
File and line references identify the implementation entry points and are
intended to be read together with the surrounding source.

CalculiX uses a shared-array architecture. The top-level driver allocates and
owns arrays for coordinates, element topology, loads, constraints, materials,
state variables, sparse matrices, and results. C and Fortran routines pass
those arrays by address across the solver pipeline. The same representation is
used by the standalone executable and, with persistent state, by
`CalculiXstep.c`.

## 1. End-to-end algorithm

The normal execution path is:

```text
job.inp
  -> readinput
  -> initialize counters and flags
  -> allocate shared arrays
  -> interpret *STEP and analysis keywords
  -> build constraints/contact structures
  -> create sparse equation structure
  -> assemble element matrices and right-hand sides
  -> solve or iterate the equation system
  -> update displacements, temperatures, stresses, and state variables
  -> recover and write .dat/.frd/restart results
```

The C entry point in `src/CalculiX.c` parses `-i jobname`, calls `readinput`,
initializes with `ini_cal`, allocates through the Fortran allocation routine,
and repeatedly enters the step interpreter. The analysis dispatch later calls
linear static, nonlinear, dynamic, modal, coupled-physics, or specialized
drivers (`src/CalculiX.c:35`, `src/CalculiX.c` around the analysis dispatch at
lines 1264-1296).

## 2. Input normalization and step interpretation

### 2.1 Input deck storage

`readinput` does more than read lines. It creates an in-memory representation
of the input deck and indexes keyword positions so later routines can process
cards in the required order (`src/readinput.c:25`).

The implementation:

1. Allocates a character buffer for the deck, an integer pointer table, and
   an integer token/index array.
2. Opens `<jobname>.inp`.
3. Reads lines, handling nested `*INCLUDE` files.
4. Removes blanks and line endings for ordinary keyword/data lines.
5. Preserves heading text differently so user headings are not normalized as
   keywords.
6. Ignores blank lines and `**` comments.
7. Converts keyword text to uppercase while preserving filename-like values
   after parameters such as `FILE=`, `INPUT=`, `OUTPUT=`, and `FILENAME=`.
8. Records keyword locations in `ipoinp` and line/token information in `inp`
   and `ipoinpc`.

This is why later Fortran routines can call `getnewline` and obtain the next
logical input card without reopening and reparsing the original file.

### 2.2 Initialization and allocation

`ini_cal` establishes the initial zero/empty state for counts and flags:
loads, boundary conditions, MPCs, materials, nodes, elements, contact,
thermal/network data, sensitivity options, and output state
(`src/ini_cal.c:22`).

`calinput` is the central step interpreter. Its header documents the analysis
method codes, perturbation/nonlinearity flags, step number, and status
conventions (`src/calinput.f:19-100`). It maintains flags such as whether
boundary conditions, loads, output requests, contact definitions, and
couplings have already occurred. Those flags control `OP=NEW` behavior and
keyword ordering.

The important algorithmic consequence is that the input deck is not converted
into a high-level object graph. Instead, keyword handlers populate shared
arrays and counters. For example:

- `materials.f` allocates a material slot and stores its name before the
  material property cards are processed (`src/materials.f:19-70`).
- `matrixassembles.f` creates a user-element-like topology and stores stiffness
  and optional mass filenames as material-associated data
  (`src/matrixassembles.f:19-150`).
- `statics.f` and `dynamics.f` parse time controls and solver selection, then
  set `nmethod`, `isolver`, and increment limits
  (`src/statics.f:19-120`, `src/dynamics.f:19-130`).

### 2.3 Static and dynamic keyword parameters

For `*STATIC`, the parser recognizes the solver, direct-solver flag, time
reset, and total-time options. It reads initial increment, step time, minimum
increment, maximum increment, and a CFD increment. It validates negative and
out-of-range values and supplies defaults for nonlinear cases
(`src/statics.f:65-270`).

For `*DYNAMIC`, the parser defaults to implicit dynamics, sets the Hilber-Hughes-
Taylor-like parameters `alpha(1)` and `alpha(2)`, and permits `EXPLICIT`.
It rejects massless contact for implicit dynamics, selects the requested sparse
solver, and applies separate increment rules for implicit and explicit modes
(`src/dynamics.f:55-275`).

These routines configure the numerical driver; they do not themselves perform
the full time/increment solve.

## 3. Finite-element integration and element matrices

### 3.1 Element dispatch

`mafillsm` is the Fortran element-loop routine for filling the sparse stiffness
matrix and related vectors (`src/mafillsm.f:19`). Its loop:

1. Iterates over a range of elements.
2. Skips deleted elements and fluid-only entries when appropriate.
3. Determines the element family and node count from `lakon`.
4. Selects degrees of freedom and special handling for solid, spring, mass,
   user, contact, and mortar elements.
5. Computes body-force data where needed.
6. Calls the element routine, normally `e_c3d` for 3-D solids or a specialized
   `e_*` routine for another formulation.
7. Maps local matrix/vector entries into global equation positions.
8. Applies MPC and boundary-condition transformations.
9. Accumulates contributions into the sparse arrays.

The C wrapper `mafillsmmain` determines thread count, partitions element
ranges, allocates per-thread matrix/vector work arrays, launches worker
threads, and reduces the thread-local contributions
(`src/mafillsmmain.c:48-220`). It then calls Fortran `mafillsm`
(`src/mafillsmmain.c:438`) and performs the force assembly path when required.

### 3.2 Shape functions, Jacobians, and quadrature

`e_c3d` is the representative 3-D element kernel
(`src/e_c3d.f:19`). Its header identifies the matrix modes:

- `nmethod=0`: positive-Jacobian check.
- `nmethod=1`: stiffness matrix and right-hand side.
- `nmethod=2`: stiffness and mass matrix.
- `nmethod=3`: static and buckling stiffness.
- `nmethod=4`: stiffness and mass matrix for dynamic use.

The routine determines the topology from `lakonl`, chooses the appropriate
integration scheme, and calls element-specific shape routines such as
`shape8h`, `shape20h`, `shape10tet`, `shape4tet`, `shape15w`, and `shape6w`
(`src/e_c3d.f:499-700`).

At each integration point the algorithm conceptually performs:

```text
reference coordinates (xi, eta, zeta)
  -> shape functions N and derivatives dN/d(xi,eta,zeta)
  -> physical Jacobian J
  -> determinant/integration weight
  -> strain-displacement operator B
  -> constitutive response sigma, tangent D
  -> local stiffness B^T D B det(J) w
  -> local internal force B^T sigma det(J) w
```

The exact implementation differs by element type, geometry, thermal mode,
composite layer, and analysis method. `e_c3d` stores local matrices in `s` and
mass-related data in `sm`; it also handles body forces, thermal terms,
prestress/buckling terms, contact-face terms, and output state.

The Jacobian check is not merely diagnostic. A nonpositive determinant causes
the surrounding matrix-fill logic to mark the analysis method as invalid,
preventing the solver from proceeding with an inverted element
(`src/mafillsmmain.c` allocation and `nmethod` handling around lines 180-220;
`src/e_c3d.f` element integration around lines 640-700).

### 3.3 Local-to-global assembly

For an element with local degrees of freedom `a` and `b`, the assembled
matrix contribution is conceptually:

```text
K[I(a), I(b)] += Ke[a,b]
fint[I(a)]     += fe[a]
```

The actual mapping is mediated by `nactdof`, `icol`, `jq`, and `irow`. These
arrays encode active equation numbers and the sparse column structure rather
than storing a dense global matrix.

`mafillsm` calls `mafillsmmatrix` when MPC-related matrix modifications are
needed (`src/mafillsm.f:304-325`). Boundary conditions and MPCs therefore
affect both the right-hand side and the matrix representation before the
external solver sees it.

## 4. Constitutive update and state variables

### 4.1 Material data versus material response

`materials.f` only registers material identities and enforces placement before
step definitions. Material property cards populate the arrays later consumed
by the element and constitutive routines. The material name in `matname` is
used as a dispatch key.

`e_c3d` selects the material associated with an element through `ielmat`,
obtains orientation through `ielorien`/`orab`, and distinguishes hyperelastic
or layered cases before evaluating integration-point behavior
(`src/e_c3d.f:280-430`).

### 4.2 Integration-point state

The principal constitutive state arrays are:

- `xstateini`: state at the beginning of the increment/step.
- `xstate`: current integration-point state.
- `sti`/`stx`: stress-like integration-point storage.
- `xstiff`: cached or derived stiffness information in selected paths.
- `ener` and related arrays: energy quantities.

The element routine passes strain-like vectors, temperature/time values,
orientation data, state arrays, and material property arrays to the material
dispatcher.

### 4.3 Material dispatcher

`umat_main` dispatches by the material-name prefix and forwards the same
constitutive interface to built-in or compatibility material routines
(`src/umat_main.f:19-43` and branches through approximately line 240).
Supported branches include ABAQUS-style models, creep, hyperelastic/ideal-gas
models, tension/compression-only models, crystal models, linear elastic models,
and `USER`.

The common response contract is visible in the argument list:

```text
strain/deformation inputs
  + old state and time increment
  + material constants and orientation
  -> stress vector `stre`
  -> tangent/stiffness `stiff`
  -> updated state `xstate`
  -> optional suggested increment `pnewdt`
```

The tangent is essential to nonlinear assembly because it supplies the local
linearization used in the global Newton correction. The state arrays are
updated per integration point, not once per element.

## 5. Sparse matrix construction and equation solving

### 5.1 Sparse representation

The main mechanical sparse representation separates diagonal values (`ad`) from
off-diagonal values (`au`). `icol`, `jq`, and `irow` describe the active
equation graph and nonzero locations. Additional arrays (`adb`, `aub`) hold a
mass or secondary matrix when an analysis requires one.

`mafillsmmain` first determines the number of active nonzeros and allocates
matrix storage. It can allocate separate stiffness, mass, external-force, and
dynamic work arrays depending on `mass`, `stiffness`, `buckling`, and `rhsi`
flags (`src/mafillsmmain.c:150-220`).

The matrix is assembled in the CalculiX sparse format, then converted by the
selected backend adapter. The representation is therefore an internal
contract: solver adapters must interpret the same `icol/irow/jq/ad/au`
semantics or explicitly convert them.

### 5.2 Iterative lumped solve

`solveeq` implements a simple iterative correction against a lumped diagonal
system (`src/solveeq.c:24`). It first computes:

```text
sol[i] = b[i] * adl[i]
```

where `adl` is the lumped inverse/scaling. For each additional iteration it
computes the sparse matrix-vector product through Fortran `spmv`, then applies:

```text
sol[i] = (b[i] - (A sol)[i]) * adl[i]
```

This is a fixed-point/Jacobi-like iteration. It is not the same path as the
direct sparse solvers and should be understood as a lightweight iterative
option, especially in coupled fluid/network routines that call it directly.

### 5.3 SPOOLES adapter

`spooles_factor` converts the CalculiX matrix into a SPOOLES `InpMtx`
representation (`src/spooles.c:459`). The conversion supports several input
formats:

- sparse lower triangular storage,
- sparse full storage,
- dense storage,
- sparse full storage with row/column conventions,
- row-wise CFD storage.

For a generalized shifted system, the adapter inserts `A - sigma B` by
subtracting the secondary matrix values while constructing the SPOOLES input.
It then changes storage mode and calls either the serial or multithreaded
factorization path.

`spooles_solve` wraps the right-hand side in a SPOOLES dense matrix, solves
using the factorization retained in `pfi`, copies the result back into the
CalculiX vector, and frees the temporary dense result
(`src/spooles.c:729-782`).

The PARDISO and PaStiX adapters follow the same architectural contract but
delegate conversion/factorization/solve details to their respective libraries
(`src/pardiso.c:37`, `src/pastix.c:179`).

## 6. Linear static analysis

`linstatic.c` is the direct linear-static driver. Its core pipeline is:

1. Allocate sparse matrix and external-force arrays.
2. Call `mafillsmmain` to assemble stiffness and load contributions.
3. If the assembled problem is asymmetric, expand/rebuild the appropriate
   matrix representation through `mafillsmasmain`.
4. Construct the right-hand side as external force minus internal force:

   ```text
   b = fext - f
   ```

5. Factor the matrix with SPOOLES, PARDISO, or PaStiX according to
   `isolver`.
6. Solve for the active degrees of freedom.
7. Call `results` to recover displacements, stresses, forces, energies, and
   requested output.

The assembly and right-hand-side construction are visible in
`src/linstatic.c:380-470`; backend selection is in the subsequent factorization
and solve branches around `src/linstatic.c:480-730`.

The linear path factors once because the matrix does not change during the
single linear solve. Repeated load cases or substructure operations may reuse
that factorization with different right-hand sides.

## 7. Nonlinear static analysis

### 7.1 Newton correction structure

`nonlingeo.c` is the nonlinear mechanical driver. Its repeated solve cycle
uses the current state to:

1. Update active constraints/MPCs and, where required, contact structures.
2. Recover stresses and internal forces from the current displacement estimate
   through `results`.
3. Assemble the tangent matrix and load vectors through `mafillsmmain`.
4. Form a residual/right-hand side from external and internal forces.
5. Factor and solve the tangent system with the selected backend.
6. Add the correction to the current solution.
7. Recompute results and convergence measures.
8. Accept the increment, reduce the increment, or stop on convergence/error.

Representative calls to `results` and `mafillsmmain` occur at
`src/nonlingeo.c:990-1090` and again in later increment/iteration paths.
Solver selection and factorization branches occur around
`src/nonlingeo.c:1260-1500` and `src/nonlingeo.c:3030-3160`.

In mathematical form, the intended local linearization is:

```text
R(u) = Fext(u) - Fint(u)
K_t(u_k) * Delta-u = R(u_k)
u_(k+1) = u_k + Delta-u
```

The source does not expose this as a small abstract Newton class. Instead,
`vold`, `v`, `f`, `fext`, `stx`, `xstate`, and the sparse arrays carry the
iteration state through large calls.

### 7.2 Increment control

The nonlinear driver works in increments. `dyna.c` provides a visible example
of the shared increment-control pattern: it advances `iinc`, computes `time`,
`dtime`, and relative step time, updates loads, and calls `checkinclength` when
automatic increment control is enabled (`src/dyna.c:1260-1340`).

For static nonlinear calculations, `*STATIC` supplies initial/minimum/maximum
increment sizes. A failed or difficult iteration can force a smaller increment;
an accepted increment may permit growth subject to the configured limits.
The exact acceptance policy is implemented in the driver and its convergence
helpers rather than in `statics.f`, which primarily parses the controls.

## 8. Contact algorithms

### 8.1 Penalty/contact-element path

`contact.c` selects the contact construction path based on `mortar`:

- `mortar <= 0`: update master-surface normals/bounding planes and call
  `gencontelem_n2f` for node-to-face contact.
- `mortar == 1`: call `gencontelem_f2f` for face-to-face mortar contact.

Before generation, it sizes temporary arrays from the largest tie/contact
triangle list. After generation, it frees all temporary geometry and index
arrays (`src/contact.c:23-105`).

The node-to-face path prepares common master-triangle planes. The generated
contact elements are then treated by the normal element/matrix assembly
pipeline, allowing contact stiffness and forces to enter the global system.

### 8.2 Dual mortar path

`contactmortar.c` implements the dual-mortar transformation and explicitly
references Algorithm 2 of the cited thesis in its source comment
(`src/contactmortar.c:23-36`).

Its source-visible stages are:

1. Build `islavactdof`, coupling active equation degrees of freedom to slave
   and master nodes.
2. Once per increment, update master-triangle geometry and bounding/search
   data.
3. On the first increment/iteration, generate an initial active set analogous
   to node-to-surface penalty contact.
4. Count active, inactive, no-gap, sticking, and non-mortar nodes.
5. Compute slave normals and tangents with `nortanslav`.
6. Locate matched slave/master integration points and construct the coupling
   data used by the mortar matrices.
7. Build or update coupling matrices for slave/master displacement, contact
   constraints, and traction terms.
8. Reduce or transform the global system using those matrices.

The source states that normals, tangents, and surface segmentation are updated
once per increment while a combined fixed-point/Newton approach handles the
remaining iterations (`src/contactmortar.c:100-160`). This is an important
performance and algorithmic distinction: geometry/search data is not rebuilt
at every Newton iteration.

The mortar routine carries sparse blocks named `B`, `D`, transformed variants,
and inverses. These represent coupling and constraint operators rather than
the ordinary element stiffness matrix. Their dimensions and construction are
visible in the long argument list and subsequent allocations in
`src/contactmortar.c`.

## 9. Dynamic time integration

### 9.1 Driver setup

`dynamics.f` parses the dynamic procedure and chooses implicit or explicit
mode. `dyna.c` performs the actual time stepping. Its major loop continues
while the normalized step time has not reached the end
(`src/dyna.c:1260` onward).

At each increment the driver:

1. Saves starting time and, when requested, energy/stress state.
2. Increments the step counter and checks the maximum-increment limit.
3. Adjusts the increment if automatic control is active.
4. Computes `time`, `dtime`, and relative time.
5. Evaluates load differences using `temploaddiff`.
6. Builds the instantaneous RHS with `rhsmain`.
7. Applies corrections for nonzero prescribed boundaries.
8. Uses the mass/stiffness system to obtain the new response.
9. Calls `results` to update stresses, forces, and output.
10. Advances state and repeats.

### 9.2 Implicit dynamics

For implicit dynamics, the code combines mass and stiffness terms into an
effective system. The source comments and matrix operations in
`src/nonlingeo.c:2520-2585` show the form:

```text
LHS = M/dt + (alpha*M + beta*K)/2
RHS = M/dt - (alpha*M + beta*K)/2
```

The actual arrays are `ad/a u` and `adb/aub` for stiffness and mass-like
components. Contact degrees of freedom may be reduced before factorization.
The effective LHS is then factored once for the operation and used for the
increment solve.

The dynamic parser constrains `alpha(1)` to the range `[-1/3, 0]` and uses
`alpha(2)` for a relative-to-absolute option
(`src/dynamics.f:80-125`). These parameters control numerical damping and
time integration behavior in the driver.

### 9.3 Explicit dynamics

For explicit dynamics, the source does not solve a full effective stiffness
system at each step. It computes a stable time increment from volumetric
stability/Courant considerations using `calcstabletimeincvol`
(`src/nonlingeo.c:2520-2560`), with optional mass-scaling modes.

With no selective mass scaling, the mechanical acceleration-like update is
formed componentwise from the force imbalance and lumped mass:

```text
b[i] = (fext[i] - f[i]) / mass_diagonal[i]
```

This is visible in `src/nonlingeo.c:2860-2890`. With selective mass scaling,
the code instead factors/solves the mass system through the selected backend
before continuing.

## 10. Eigenvalue and modal analysis

`arpack.c` wraps ARPACK's reverse-communication iteration
(`src/arpack.c:44`). It sets up a generalized eigenproblem using the assembled
stiffness and mass operators.

### 10.1 Problem setup

The routine obtains the requested number of modes from `mei`, the Arnoldi
subspace size from `mei`, and tolerance from `fei` (`src/arpack.c:183-186`).
It assembles the mechanical matrices through `mafillsmmain` and factors the
operator needed by the selected spectral transformation
(`src/arpack.c:477-529`, `src/arpack.c:620-683`).

For a generalized problem, the conceptual target is:

```text
K phi = lambda M phi
```

The source uses ARPACK's reverse-communication interface, so ARPACK asks the
caller to perform matrix-vector or transformed solves rather than receiving
explicit matrix callbacks.

### 10.2 Reverse communication

The routine initializes ARPACK parameters, allocates the Krylov basis, and
loops while `ido` is `-1`, `1`, or `2` (`src/arpack.c:694-749`).

Within that loop:

- `ido == -1` requests the initial operator action.
- `ido == 1` requests the main transformed operator action, usually involving
  the previously factored matrix.
- `ido == 2` requests a mass-matrix or secondary operator action.

The caller performs the requested operation on the work arrays and calls
`dnaupd` or `dsaupd` again (`src/arpack.c:749-850`). Once the iteration reports
convergence, `dneupd` or `dseupd` extracts eigenvalues and eigenvectors
(`src/arpack.c:855-883`).

The routine then writes eigenvalues, computes modal quantities, and calls
`results` for mode-specific output. Thus modal output is not a separate
postprocessor; it reuses the ordinary result-recovery path.

## 11. Result recovery

`results.c` converts the solved primary variables into physical output
(`src/results.c:46`). Its documented modes distinguish whether the primary
field is already known, whether it should be computed from the system
solution, and whether output is requested.

The high-level result algorithm is:

1. `resultsini` examines requested labels and determines which derived
   quantities are needed: internal forces, external forces, heat quantities,
   Cauchy stresses, integration-point variables, and output fields
   (`src/results.c:170`).
2. The element range is partitioned across worker threads.
3. `resultsmech` or `resultstherm` evaluates element-level stresses, strains,
   heat fluxes, material tangents, and nodal forces
   (`src/results.c:180-250`, `src/results.c:493-525`).
4. Thread-local internal-force arrays are reduced into the global `fn`.
5. Requested nodal and integration-point fields are sent to `resultsprint`.
6. `frd` and related writers serialize visualization data.

This separation matters for debugging: a correct linear solve can still expose
a result regression if `resultsmech`, `resultsprint`, `frd`, or an output-label
selection path changes.

## 12. Parallel decomposition

The core parallel pattern is element-range partitioning:

1. Determine the available CPU count from `NUMBER_OF_CPUS`,
   `CCX_NPROC_*`, or `OMP_NUM_THREADS`.
2. Limit the requested count to available CPUs and, where applicable, the
   number of elements.
3. Split the element interval into per-thread ranges.
4. Allocate thread-local arrays for quantities that would otherwise race.
5. Run element assembly or result recovery in parallel.
6. Join threads and reduce thread-local arrays into the shared result.

`mafillsmmain` implements this for stiffness/mass assembly
(`src/mafillsmmain.c:100-220`). `results` implements the same pattern for
stress and internal-force recovery (`src/results.c:100-220`).

The design avoids concurrent writes to the main sparse arrays by assembling
into per-thread buffers and reducing them. Any new parallel element operation
must preserve that ownership rule.

## 13. Algorithm-to-source index

| Algorithm | Main implementation | Supporting sources |
| --- | --- | --- |
| Input normalization | `src/readinput.c` | `src/calinput.f`, `src/inputerror.f` |
| Step/analysis selection | `src/calinput.f` | `src/statics.f`, `src/dynamics.f` |
| Element integration | `src/e_c3d.f` | `src/shape*.f`, `src/e_*.f` |
| Constitutive dispatch | `src/umat_main.f` | `src/material*.f`, `src/umat*.f` |
| Sparse assembly | `src/mafillsmmain.c`, `src/mafillsm.f` | `src/matrixassembles.f`, `src/matrixstorage.c` |
| Linear static | `src/linstatic.c` | solver adapters, `src/results.c` |
| Nonlinear static | `src/nonlingeo.c` | `src/contact*.c/f`, `src/results.c` |
| Contact generation | `src/contact.c` | `gencontelem_n2f`, `gencontelem_f2f` |
| Dual mortar | `src/contactmortar.c` | mortar sparse block routines |
| Dynamic stepping | `src/dyna.c` | `src/dynamics.f`, `src/nonlingeo.c` |
| Eigenvalue iteration | `src/arpack.c` | `src/arpackbu.c`, `src/arpackcs.c` |
| Direct sparse solve | `src/spooles.c`, `src/pardiso.c`, `src/pastix.c` | `CMakeLists.txt` solver configuration |
| Result recovery | `src/results.c` | `src/resultsmech.f`, `src/resultstherm.f`, `src/frd.c` |

## 14. Reading the algorithms safely

Several source-level details are easy to misinterpret:

- `statics.f` and `dynamics.f` primarily parse procedure parameters; the
  actual solve loops are in C drivers such as `linstatic.c`, `nonlingeo.c`,
  and `dyna.c`.
- `e_c3d` is representative, not universal. Shells, beams, fluids, network
  elements, user elements, and electromagnetic elements use specialized
  routines.
- `ad/au` and `adb/aub` do not always mean the same physical pair in every
  specialized driver. Always inspect the caller's `mass`, `stiffness`,
  `buckling`, `inputformat`, and `sigma` flags.
- Contact has multiple formulations. Node-to-face penalty contact, face-to-face
  mortar contact, tied contact, and contact elements do not share identical
  update equations.
- Result recovery may call constitutive routines again to calculate derived
  stresses and output quantities. It is not necessarily a passive dump of the
  solver vector.
- The source supports both LP64 and ILP64 builds. Index types and external
  solver ABI assumptions are part of the algorithm's executable contract.
