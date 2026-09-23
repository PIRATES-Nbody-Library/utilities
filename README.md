# Utilities

This suite contains auxiliary software to set up, organize, post-process, and
analyze the results from SWIFT N-body simulations.

## Prerequisites

- A Fortran compiler (`gfortran`, `ifort`, or `ifx`)
- A built copy of the
  [Nbody_simulation_suite](https://github.com/PIRATES-Nbody-Library/Nbody_simulation_suite),
  with `libswift.a` compiled

## Installation

### Clone from GitHub

```
git clone git@github.com:PIRATES-Nbody-Library/utilities.git
```

### Configure and generate files

From the `utilities/` directory, run `setup.sh`.  By default it looks
for the simulation suite at `../Nbody_simulation_suite/swift_pirates`:

```
./setup.sh
```

To point to a different location:

```
./setup.sh --swift-dir=/path/to/swift_pirates
```

### Build all utilities

```
./setup.sh --build
```

Or build a single utility:

```
./setup.sh --build=gen_all_multi
```

### Debug build

```
./setup.sh --debug --build
```

### Clean generated files

```
./setup.sh --clean
```

## Directory Layout

After running `setup.sh`:

```
utilities/
├── setup.sh                        # Configuration script
├── TODO.md
├── README.md
├── .gitignore
├── bin/                             # Compiled executables
├── fortran/
│   ├── @make_utilities.template     # Build script template
│   ├── @make_utilities              # Generated build script
│   ├── *.f                          # Fortran source files (tracked)
│   ├── gen_all_multi.template.f     # Template (tracked, needs setup.sh)
│   ├── gen_multi_reprise.template.f # Template (tracked, needs setup.sh)
│   ├── gen_all_multi.f              # Generated from template (ignored)
│   └── gen_multi_reprise.f          # Generated from template (ignored)
├── sh/
│   ├── gen_all_multi.template.sh    # Shell script template (tracked)
│   ├── gen_all_multi.sh             # Generated shell script (ignored)
│   ├── paramsmulti.template.in      # Input file template (tracked)
│   └── paramsmulti.in               # Generated input file (ignored)
└── workspace/                       # Runtime directories (ignored)
    ├── working/
    ├── TMP/
    └── multi/
```

## Options

Run `./setup.sh --help` for all options including SBATCH parameters
for HPC job submission and custom directory paths.

## Status

This repository is in early development.
