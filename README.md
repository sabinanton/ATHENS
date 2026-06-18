# ATHENS

Aerodynamic Tracer for High-altitude Environment Numerical Simulations (ATHENS) is an open-source C++17 ray-tracing solver for aerodynamic modelling of Resident Space Objects (RSOs) in free-molecular, hyperthermal flow. It is intended for orbital-aerodynamics cases where intermolecular gas collisions can be neglected on spacecraft length scales, allowing the flow problem to be treated through particle ray tracing rather than full DSMC or TPMC simulation.

ATHENS imports spacecraft geometries from STL files, places them inside a cuboidal virtual wind-tunnel domain, and spawns gas rays from the domain boundaries according to the prescribed bulk and thermal flow conditions. Each ray represents a molecular mass flux with a velocity composed of a free-stream component and a Maxwellian thermal component. Ray-surface intersections are detected on triangular surface elements, after which reflected rays are generated using a selected gas-surface interaction scattering kernel.

The code currently supports the Diffuse Re-emission with Incomplete Accommodation (DRIA) and Cercignani-Lampis-Lord (CLL) scattering kernels. These models are used to compute the momentum transferred from incident and reflected particles to the spacecraft surface, including repeated reflections for concave geometries.

ATHENS computes global aerodynamic force and moment coefficients, including drag, lift, and side-slip components. It can also estimate surface-level pressure and shear coefficient distributions on the triangular elements of the imported geometry.

The code can be run locally or on an HPC cluster using MPI.

## Repository structure

```text
.
├── athens              # Compiled executable, generated after compilation
├── src/                # Source files
├── code_files/         # Additional code and configuration files
├── test_project/       # Example project directory
├── compile.sh          # Compilation script
├── run.sh              # MPI run script
├── postprocess.sh      # Post-processing script
├── plot.sh             # Plotting script
├── clean.sh            # Cleanup script
├── LICENSE
└── README.md
```

## Requirements

ATHENS requires:

* A C++ compiler
* MPI
* `mpirun`
* Bash

On an HPC system, load the required compiler and MPI modules before compiling or running the code.

Example:

```bash
module load gcc
module load openmpi
```

The exact module names depend on the system.

## Compilation

Compile the code from the repository root:

```bash
./compile.sh
```

If compilation is successful, the executable is created as:

```text
./athens
```

## Running ATHENS

ATHENS is run using:

```bash
./run.sh -project=<project_name> -np=<number_of_processes> -id=<simulation_id>
```

Example:

```bash
./run.sh -project=test_project -np=16 -id=0
```

This runs ATHENS in parallel using MPI:

```bash
mpirun -np <number_of_processes> ./athens --project=<project_name> --id=<simulation_id>
```

## Batch mode

To run in batch mode, add the `-batch` flag:

```bash
./run.sh -project=<project_name> -np=<number_of_processes> -id=<simulation_id> -batch
```

Example:

```bash
./run.sh -project=test_project -np=16 -id=0 -batch
```

This passes the `--batch` flag to the ATHENS executable.

## Command-line arguments

The `run.sh` script uses the following arguments:

| Argument                    | Description                            |
| --------------------------- | -------------------------------------- |
| `-project=<project_name>`   | Name of the project directory          |
| `-np=<number_of_processes>` | Number of MPI processes                |
| `-id=<simulation_id>`       | ID of the simulation run               |
| `-batch`                    | Optional flag for batch-mode execution |

Internally, these are passed to the ATHENS executable as:

| ATHENS argument            | Description                   |
| -------------------------- | ----------------------------- |
| `--project=<project_name>` | Name of the project directory |
| `--id=<simulation_id>`     | ID of the simulation run      |
| `--batch`                  | Enables batch-mode execution  |

## Project directories

Each simulation is associated with a project directory. The project directory contains the input files and stores the output files for a given set of simulations.

A typical run is:

```bash
./run.sh -project=test_project -np=16 -id=0
```

Here:

* `test_project` is the project directory
* `16` is the number of MPI processes
* `0` is the simulation ID

## Output

Simulation output is written inside the selected project directory.

For batch simulations, ATHENS may create a batch-results directory such as:

```text
<project_name>/Batch_Results/
```

Post-processing and plotting can then be performed with:

```bash
./postprocess.sh
./plot.sh
```

Check the relevant scripts before running them, since they may assume specific project names, file paths, or output formats.

## Gas species

The code includes constants for common thermospheric gas species:

| Species            | Symbol |
| ------------------ | ------ |
| Atomic oxygen      | AO     |
| Molecular oxygen   | O2     |
| Helium             | He     |
| Molecular nitrogen | N2     |
| Argon              | Ar     |

These constants are used together with the gas constant, Avogadro constant, and Boltzmann constant in the gas-surface interaction calculations.

## Running on DelftBlue

To run on DelftBlue, first connect to the login node:

```bash
ssh <netid>@login.delftblue.tudelft.nl
```

Move to the ATHENS repository:

```bash
cd athens_4tu
```

Compile:

```bash
./compile.sh
```

Run a simulation using MPI:

```bash
./run.sh -project=<project_name> -np=<number_of_processes> -id=<simulation_id>
```

Run a batch simulation:

```bash
./run.sh -project=<project_name> -np=<number_of_processes> -id=<simulation_id> -batch
```

Use realistic wall-time requests when submitting jobs on the cluster.

## Cleaning generated files

To remove generated files, use:

```bash
./clean.sh
```

Check the contents of `clean.sh` before running it, since cleanup scripts may remove compiled files or simulation output.

## Typical workflow

```bash
# 1. Compile
./compile.sh

# 2. Run a test case in parallel
./run.sh -project=test_project -np=16 -id=0

# 3. Run a batch case
./run.sh -project=test_project -np=16 -id=0 -batch

# 4. Post-process results
./postprocess.sh

# 5. Plot results
./plot.sh
```

## Notes

* Run all scripts from the repository root unless stated otherwise.
* Ensure that `athens` exists before running `run.sh`.
* Ensure that `mpirun` is available in the current environment.
* On HPC systems, compile and run using compatible compiler and MPI modules.
* Keep large simulation output on scratch storage when running on a cluster.
* Do not assume scratch storage is backed up.
