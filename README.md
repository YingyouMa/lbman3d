# lbman3d
Lattice Boltzmann Method-based solver for 3D Active Nematics.

The flow equation is solved using a D3Q15 scheme. The Q-tensor equation is solved using a finite-difference scheme.

## Dependencies

| Dependency | Version | Notes |
|---|---|---|
| CMake | ≥ 3.23 | Build system |
| C++ compiler | C++23 | GCC 13+ or Clang 17+ recommended |
| OpenMP | — | Usually bundled with the compiler |
| HDF5 | any recent | C library only |
| [kokkos/mdspan](https://github.com/kokkos/mdspan) | `stable` | Fetched automatically by CMake |

## Solver Workflow

This repository now follows a template-case workflow similar to the HPCC multigrid solver.

### Workspace layout

- `src/`
  Shared solver source code.
- `template/`
  Master case template. Copy this to start a new run.
- `project_runs/`
  Actual simulation cases copied from `template/`.
- `readme/`
  Workflow notes.

### Create a new case

From the repository root:

```bash
mkdir -p project_runs
cp -r template project_runs/<case_name>
cd project_runs/<case_name>
```

For a normal run, edit only these files inside the copied case:

- `boundary_config.h`
- `init_config.h`
- `parameters.h`
- `run_config.h`

Everything under `src/` is shared solver code and should not need changes for routine simulations.

## Build And Run

### Local case build

From inside a case directory:

```bash
cmake -S . -B build-case
cmake --build build-case -j$(nproc)
./build-case/main
```

The case-local `CMakeLists.txt` points back to the shared solver source and injects the case configuration headers automatically.

### HPCC batch run

From inside a case directory:

```bash
./submit.sh
```

This submits `example.sh`, which configures, builds, and runs the solver from that case directory.

## Output

Outputs are written relative to each case directory:

- `data/`
- `restart/`
- `lbm.log`
- `parameters.auto.md`
- `output.resume`
- `errors.resume`
- `slurm-<jobid>.out`
- `slurm-<jobid>.err`

## Notes

The repository root can still be built directly for development, using the default case configuration from `src/default_case_config.h`.
