# NPDECODES

Reviewing the problems collection has 3 stages:
1. Polishing the code
2. Verifying consistency with the pdf document (creating a consistency txt file)
3. Modifying the latex associated to the pdf document to fix the inconsistencies and
replace the old links by new ones pointing to the correct /NPDECODES/homeworks folder on github

- A white check mark :white_check_mark: indicates that a stage is in progress.
- A green check mark :heavy_check_mark: indicates that a stage is completed.
- Under 'Assignee for current stage' is found the name of the assistant currently working on the problem.

| # | Problem name | Polished | Verified | Latex | Assignee (stage in progress)|
| --- | --- | --- | --- | --- | --- |
| 2-2 | `TransformationOfGalerkinMatrices` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |
| 2-4 | `LinearFE1D` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |
| 2-5 | `SimpleLinearFiniteElements` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |
| 2-6 | `IncidenceMatrices` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |
| 2-7 | `LengthOfBoundary` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |
| 2-8 | `ElementMatrixComputation` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | |
| 2-9 | `LFPPDofHandling` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | |
| 2-10 | `ProjectionOntoGradients` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | |
| 2-12 | `TestQuadratureRules` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |
| 2-13 | `ParametricElementMatrices` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |
| 2-14 | `NonConformingCrouzeixRaviartFiniteElements` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |
| 2-15 | `RegularizedNeumannProblem` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | |
| 3-1 | `AvgValBoundary` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | |
| 3-2 | `DebuggingFEM` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |
| 3-3 | `PointEvaluationRhs` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | |
| 3-4 | `UnstableBVP` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |
| 3-5 | `ErrorEstimatesForTraces` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | |
| 3-7 | `MaximumPrinciple` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | |
| 3-8 | `OutputImpedanceBVP` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark:  |  |
| 6-1 | `RadauThreeTimestepping` | :heavy_check_mark: | :heavy_check_mark: |  |  |
| 6-2 | `SDIRKMethodOfLines` | :heavy_check_mark: | :heavy_check_mark: |  |  |
| 6-4 | `1DWaveAbsorbingBC` | :heavy_check_mark: |:heavy_check_mark: |  |  |
| 6-5 | `SymplecticTimesteppingWaves` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | |
| 6-6 | `BoundaryWave` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |
| 8-1 | `BurgersEquation` | :heavy_check_mark: | :heavy_check_mark: |  |   |
| 8-2 | `EngquistOsherNumericalFlux` | :heavy_check_mark: | :heavy_check_mark: |  |  |
| 8-3 | `FiniteVolumeSineConsLaw` | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |
| 8-6 | `CLEmpiricFLux` | :heavy_check_mark: | :heavy_check_mark: |  |  |
| ? | `WaveABC2d` | :heavy_check_mark: |  |  |  |
| ? | `ZienkiewiczZhuEstimator` | :heavy_check_mark: |  |  |  |
| ? | `ElectrostaticForce` | :heavy_check_mark: |  |  |  |
| ? | `ExtendedMUSCL` |  |  |  |  |
| ? | `FiniteVolumeRobin` |  |  |  |  |
| ? | `IPDGFEM` | |  |  |  |
| ? | `LinFeReactDiff` | |  |  |  |

## General Remarks

* Below, the directory `.` refers to the top level directory of the repository.
* Developers only work in `./developers/<ProblemName>/`. From this, a corresponding directory `./homeworks/<ProblemName>/` for the students can be created automatically using `./scripts/deploy_npde.py`.
* Not only the solutions, but also the corresponding templates need to compile and, if possible, run without crash. So be careful when setting the solution/template tags in `developers/mastersolution/`.
* The bullets below are only a selection. If you spot additional issues, e.g. ugly or too complicated code, fix it.
* Names of .cc and .h files: For example the files in the folder `./developers/MyHomeworkProblem/mastersolution/` should be called `myhomeworkproblem_main.cc`, `myhomeworkproblem.cc`, `myhomeworkproblem.h`, `myclass.h` (lowercase, no underlines except for the main file).
* Only in Lehrfem exercises: Use `nostd::span` (C++20) instead of `ForwardIteraters` for iterating over objects contiguous in memory (used e.g. in `Mesh::Entities()`, `SubEntities()`, `DofHandler`).

## Polishing

* Add comments if necessary.
* Use `#include <cmath>` instead of `#include <math.h>`.
* Only use the `#includes "..."` you need, e.g. for `Eigen::MatrixXd`, `#include <Eigen/Core>` is sufficient (no `#include <Eigen/Dense>` needed). Likewise for the `dependencies.cmake`. However, every file (`.h` and `.cc`) should contain all the `#includes` it needs, not relying on other headers to provide them silently.
* Use forward declarations to avoid `#includes` whenever possible.
* Make the code follow the Google style (see https://google.github.io/styleguide/cppguide.html). In particular, respect the include order (https://google.github.io/styleguide/cppguide.html#Names_and_Order_of_Includes).
* Never store data in the remote repository that is generated by the code anyway (e.g. figures). The repo should only contain source files (and maybe some mesh files).
* Provide unit tests for all functions that the students have to change.
* Plotting is done as follows: The C++ code only stores the data in a `.csv` file (using `std::ofstream`). The plotting is then done by a python script that reads this `.csv` file. The python script is called within the C++ code via the `std::system(...)` function. See e.g. `BurgersEquation`.

## Homework Problems

### Setting Up a Homework Problem

Only work in `./developers/<ProblemName>/`. This folder needs to contain at least:
* `CMakeLists.txt` (same for every problem, so just copy-paste it from another problem)
* `mastersolution/` (contains you code, with tags as described below)
* `mastersolution/dependencies.cmake` (contains the info needed by the build system, just mimic the ones from other problems)
Every problem should contain unit test in the directory `mastersolution/test/` containing the following two files:
* `problemname_test.cc` (contains the unit tests for your code)
* `dependencies.cmake` (contains the info needed by the build system, just mimic the ones from other test directories)
In addition, if a problem needs to load meshes, we put them in the folder `./developers/<ProblemName>/meshes/`. Finally, the line `add_subdirectory(<ProblemName>)` has to be added in `./developers/CMakeLists.txt`.

### Solution Tags

In the files of `./developers/mastersolution/` we put the following tags
```
#if SOLUTION
  <only in mastersolution>
#else
  <only in template>
#endif
```
to indicate what belongs to mastersolution and/or template. Based on these tags, the file `./scripts/deploy_npde.py` generates a directory `./homeworks/<ProblemName>/` containg the directories `mastersolution`, `mysolution`, `temaplates` with the corresponding content. The students work exclusively in `./homeworks/<ProblemName>/`.

## New Problems

Problems PDF: https://www.sam.math.ethz.ch/~grsam/NUMPDE/HOMEWORK/NPDEProblems.pdf

* Problem 5.6: Parametric Finite Elements
* Problem 5.7: Stable Evaluation at a Point
* Problem 5.8: Trace Error Estimates **(done already?)**
* Problem 6.6: Non-linear Schrödinger Equation with Cubic Non-Linearity **(Oliver)**
* Problem 7.3: Upwind Quadrature
* Problem 7.4: Exponentially fitted upwind scheme
* Problem 7.5: Transport Problem
* Problem 7.6: Upwind Finite Volume Method

# TO DO LIST

## Missing unitests

RegularizedNeumannProblem: getGalerkinLSE(...) has no unit tests

### Not in NPDEFL_Problems
- BoundaryWave (done)
- ZienkiewiczZhuEstimator (done)
- WaveABC2D (done)

### Chapter 6
- RadauThreeTimestepping
- SDIRKMethodOfLines
- SymplecticTimesteppingWaves

