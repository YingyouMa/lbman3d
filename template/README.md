# lbman3d case template

Copy this directory into `project_runs/<case_name>` and edit only these files for normal runs:

- `boundary_config.h`
- `init_config.h`
- `parameters.h`
- `run_config.h`

Then build and run from inside the case directory:

```bash
cmake -S . -B build-case
cmake --build build-case -j$(nproc)
./build-case/main
```

For HPCC submission:

```bash
./submit.sh
```

To initialize both Q and velocity from an external file, switch `SelectedInitialCondition`
in `init_config.h` to `ExternalH5InputIC` and point `input_file` at an HDF5 file containing:

- `qxx`, `qxy`, `qxz`, `qyy`, `qyz`
- `ux`, `uy`, `uz`
- optionally `rho`
