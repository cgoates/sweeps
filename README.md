# Sweeps

[![build-and-run-tests-macos](https://github.com/cgoates/sweeps/actions/workflows/build_and_run_tests_macos.yaml/badge.svg?branch=master)](https://github.com/cgoates/sweeps/actions/workflows/build_and_run_tests_macos.yaml) [![build-and-run-tests-ubuntu](https://github.com/cgoates/sweeps/actions/workflows/ubuntu_build_and_run_tests.yaml/badge.svg?branch=master)](https://github.com/cgoates/sweeps/actions/workflows/ubuntu_build_and_run_tests.yaml)

This is a research code investigating techniques to create coarse hex layouts for fitting multipatch B-splines to be used in isogeometric analysis.
It uses foliations and harmonic functions to extend quad layouts to sweepable volumetric geometries.
By restriction to linear B-splines, it can also generate hex meshes.
The input to the algorithm is a tet mesh on the sweepable geometry, an indication of which portions of the boundary the sweep runs from and to (the source and target surfaces), and a quad layout on the source surface of the sweep which will be swept through the volume.

Papers
---
The following papers have been written based on the material in this codebase:

C. B. Goates and K. M. Shepherd, "Harmonic-based sweeps need not yield volumetric parameterizations," *Computer-Aided Design* (2025). DOI: [10.1016/j.cad.2025.104007](https://doi.org/10.1016/j.cad.2025.104007)

Presentations
---
The following presentations have been given based on the material in this codebase:

C. B. Goates and K. M. Shepherd, "Curvilinear hexahedral cell generation for swept trivariate splines using foliations," *12th International Conference on IsoGeometric Analysis* (2024).

C. B. Goates and K. M. Shepherd, "Hexahedral mesh generation for swept splines using foliations," *2025 SIAM International Meshing Roundtable* (2025).

C. B. Goates and K. M. Shepherd, "Curvilinear Hexahedral Cell Complexes for Swept Trivariate Splines via Foliations," *2025 SIAM Computational Science and Engineering* (2025).

C. B. Goates and K. M. Shepherd, "Enabling high-order computational flow simulations using foliation-informed hexahedral spline discretizations," *18th U.S. National Congress on Computational Mechanics* (2025).


Building
---
This project requires CMake 3.24 or newer.
Once you have cloned the repo, build as follows:

In the build directory of choice (e.g. `<git root>/build`) run
```
cmake <path/to/git/root> -DCMAKE_BUILD_TYPE=Release
make -j10
```
Dependencies will be downloaded and built automatically.
To run tests, run `ctest -j10` in the build directory after building.

Recreating Results
---
To recreate the results published in the journal article "Harmonic-based sweeps need not yield volumetric parameterizations" ([10.1016/j.cad.2025.104007](https://doi.org/10.1016/j.cad.2025.104007)), follow these steps:
Once the code is built and the tests are passing, run the following command in the build directory to generate the results for the bordered surface counterexample that is a union of two cylinders:
```
src/sweep output-laplace counter2 output-critical-points
```
This will generate two files: `test.vtu`, containing the tetrahedral mesh with the discrete harmonic function defined on it as a field called "laplace" and `critical_points.vtu`, which contains the critical points of the discrete harmonic function.

To generate the results for the closed surface counterexample, run the command
```
src/sweep output-laplace part_torus_in_sphere output-critical-points
```
which will generate the same files, overwriting them if they exist.
These `.vtu` files can be opened in ParaView to generate level sets using their Contour filter.

API
---
The sweeps codebase includes a work-in-progress Python API.
To use the python API, after building the codebase, include the following in your python script:

```
import sys
sys.path.insert(0, '<path/to/build/directory>/src/api')
import sweeps
```

Documentation for the sweeps api can be found in python using `help(sweeps)`.
An example use of the api can be seen in `scripts/test_sweep_param.py`, particularly in the `meshHookWithQuadMesh()` function.
