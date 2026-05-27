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
