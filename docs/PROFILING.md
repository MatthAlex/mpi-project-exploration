# Profiling

## MAQAO

[MAQAO](maqao.org) (Modular Assembly Quality Analyzer and Optimizer) is a performance analysis and optimization framework operating at binary level with a focus on core performance. Its main goal of is to guide application developpers along the optimization process through synthetic reports and hints.

The sample `template_mpi.json` under `maqao/` is used to instrument this tool.

### First steps and setup

- Grab the tar: `wget https://maqao.org/maqao_archive/maqao.x86_64.2.20.1.tar.xz`
- Unpack: `tar -xf maqao.x86_64.2.20.1.tar.xz`
- Export the binary path for ease of use: `export MAQ=path/to/maqao.x86_64.2.20.1/bin/maqao`

### Sample use

MAQAO is fairly complex in its use. Outlined here is one such command. Consult the documentation for more examples and uses.

First use `fpm` to compile and install the program with the appropriate flags: `-g -fno-omit-frame-pointer` and any other optimization options appropriate. Interfacing `maqao` and `fpm run` is non-trivial, while installing the binary is.

```bash
fpm @installgcc-maqao  # see fpm.rsp for more information
```

Then invoke `maqao`, pointing to its config file.

```bash
$MAQ oneview -R1 --config=maqao/template_mpi.json -xp=maqao/test -force-static-analysis  -maximum-threads-per-process=2
```

This will run the executable as given inside [`template_mpi.json`](../maqao/template_mpi.json) with the appropriate commands, and generate the output under `maqao/test/`.

After profiling is finished, if on a cluster, `rsync/scp` the `test/` dir locally and open it with a browser:

```bash
rsync -avz --partial user@login.hpc.cluster.ac.uk:/path/to/mpi-project-exploration/maqao/test .
firefox test/RESULTS/sendrecv_3D_one_html/index.html
```

### Caveats

The config is set up for MPI profiling on a single node, and OpenMP set to 1. Adjust the config accordingly.
