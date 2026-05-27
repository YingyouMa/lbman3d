# Solver Workflow

This workspace is organized around a template-case workflow similar to the multigrid solver.

## Workspace Layout

- `src/`
  Shared solver source code. Normal runs should not require edits here.
- `template/`
  Master case template. Copy this when starting a new run.
- `project_runs/`
  Case directories created from `template/`.
- `readme/`
  Workflow notes.

## Create a New Case

From the repository root:

```bash
mkdir -p project_runs
cp -r template project_runs/<case_name>
cd project_runs/<case_name>
```

Then edit only these files:

- `boundary_config.h`
- `init_config.h`
- `parameters.h`
- `run_config.h`

## Build And Run Locally

From inside the case directory:

```bash
cmake -S . -B build-case
cmake --build build-case -j$(nproc)
./build-case/main
```

## Submit On HPCC

From inside the case directory:

```bash
./submit.sh
```

The batch script will configure, build, and run the solver inside the case directory.

## Output Files

Normal outputs are written relative to the case directory:

- `data/`
- `restart/`
- `lbm.log`
- `parameters.auto.md`
- `output.resume`
- `errors.resume`
- `slurm-<jobid>.out`
- `slurm-<jobid>.err`
