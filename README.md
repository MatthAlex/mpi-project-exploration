# Fortran MPI exploration

Exploring a simple MPI cartesian grid with halo region updates, `mpi_f08`, `havaita`, some clever encapsulation, and general project structure and architecture.

- [Specifications](#specifications)
- [Getting Started](#getting-started)
- [Development Tools](#development-tools)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [Contact](#contact)
- [License](#license)
- [Appendix](#appendix)

## Specifications

Key features implemented:

- MPI Cartesian Grid
    - 3D topology with configurable dimensions and periodicity
    - Automatic decomposition via `MPI_Dims_create`
    - Neighbor discovery using `MPI_Cart_shift`
- Halo Region Updates
    - Blocking communication via `MPI_SendRecv`
    - Static transfer buffers for real and integer types
    - Interface-based type handling
- Boundary Conditions
    - Periodic (MPI-handled), Dirichlet, Neumann
    - Integrated with MPI topology for domain boundaries
    - Parameterized boundary types per face

### Future Improvements

- Modern MPI Interface
    - Migration to `mpi_f08` for improved type safety
    - Enhanced error handling and propagation
    - Fortran 2008 submodules for implementation hiding
- Advanced Communication
    - Non-blocking halo exchanges via `MPI_Isend`/`MPI_Irecv`
    - One-sided communication options for boundaries
    - Derived datatypes for halo regions
- Development Infrastructure
    - Parallel testing framework with `fpm`
    - Distributed logging system
    - Performance benchmarking suite

## Getting Started

### Prerequisites

- Linux OS or WSL2 via Windows
- GFortran or oneAPI Fortran compiler
- MPI/MPICH installation. Currently testing on OpenMPI v5.0.2

#### Optional

- VS Code installation
- Modern Fortran extension installation

### Step-by-step instructions

0. Load compiler modules (if using a compute cluster):

    ```sh
    module load gcc openmpi
    # or
    module load intel intelmpi
    ```

1. Create and activate a Python environment:

    ```sh
    python3 -m venv .venv
    source .venv/bin/activate  # Remember to activate your enviroment before runtime or development tasks.
    ```

2. Install the runtime and development packages:

    ```sh
    pip3 install .  # Install runtime dependencies from pyproject.toml
    pip3 install .[dev]  # Install development dependencies
    ```

3. Run the tests and main program:

    ```sh
    fpm test  # by default it uses gcc/gfortran to compile and run
    fpm run sendrecv_1D
    ```

## Development Tools

The template includes preconfigured development tools and settings for Modern Fortran development, with optimized configurations for:

- VS Code integration
- Language server features
- Code formatting
- Automated testing
- Package management

For detailed setup instructions and tool configurations, see [TOOLING.md](./docs/TOOLING.md).

## Documentation

Automated documentation is generated by using the [ford](https://github.com/Fortran-FOSS-Programmers/ford) package. For example config files and usage, see [here](https://forddocs.readthedocs.io/en/latest/index.html).

To generate dcoumentation for this sample project, and then view it on the browser, run:

```sh
ford ford.md
firefox docs/ford/index.html  # Alternatively, use your preferred browser
```

## Contributing

Contributions from the community are welcome. To contribute, consider opening an issue or pull request with changes and suggestions.

## Contact

For questions or suggestions, please contact me at [email](matt.alexandrakis@gmail.com) or open an issue.

## License

The project is operating under an [MIT](./LICENSE) license. You are free to use, modify, and distribute the code as needed for your project. Feel free to adapt and customize it to suit your requirements.

## Appendix

### Project Directory Structure

```sh
$ tree -Ia '__pycache__|.git|.pytest_cache|.venv|build|.gen*|ford'
.
├── app  # The main program driver(s) resides here
│   └── main.f90
├── docs
│   ├── MIGRATION.md
│   └── TOOLING.md
├── ford.md  # FORD config file
├── .fortls  # VSCode Modern Fortran config file
├── fpm.toml  # Fortran Package Manager config file
├── .fprettify.rc  # fprettify config file
├── .gitignore  # Git ignore list of files and directories
├── LICENSE
├── .pre-commit-config.yaml  # pre-commit config file
├── pyproject.toml  # config file
├── README.md  # you are here!
├── src  # All source code files are placed in here, except main driver
│   └── first_steps.f90
├── test  # All tests are placed in here
│   └── check.f90
└── .vscode  # Holds VSCode configs and runtime/debugging tasks
    ├── extensions.json  # simply populates the "Recommended" Extensions tab
    └── settings.json  # also referred to as "Workspace Settings (JSON)"
```

### References and Links

- [template](https://github.com/MatthAlex/fortran-project-template)
- [`fpm`](https://github.com/fortran-lang/fpm)
- [Modern Fortran extension](https://github.com/fortran-lang/vscode-fortran-support)
- [`fortls`](https://github.com/fortran-lang/fortls)
- [`fprettify`](https://github.com/pseewald/fprettify)
- [`pre-commit`](https://pre-commit.com/)
- [`ford`](https://github.com/Fortran-FOSS-Programmers/ford)
- [`uv`](https://github.com/astral-sh/uv)
