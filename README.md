# CoFlame
CoFlame simulates the sooting formation inside laminar coflow diffusion flames.

This repository contains the code version 1.8, with updated CMake build system and utilities.

Refer to the CoFlame paper for a full description.

```
@article{eaves2016coflame,
  title={CoFlame: A refined and validated numerical algorithm for modeling sooting laminar coflow diffusion flames},
  author={Eaves, Nick A and Zhang, Qingan and Liu, Fengshan and Guo, Hongsheng and Dworkin, Seth B and Thomson, Murray J},
  journal={Computer Physics Communications},
  volume={207},
  pages={464--477},
  year={2016},
  publisher={Elsevier}
}
```

DOI: 10.1016/j.cpc.2016.06.016

# Build instruction

Requirements:
 * CMake (version 3.10 or newer)
 * Fortran compiler
 * MPI library

```
# configure

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# build
cmake --build build
```

# Running examples

The `examples` directory contains a methane flame simulation `5atm`. Refer to the `run.sh` script for the steps to execute the simulation.

1. The chem

