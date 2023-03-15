# <img src=".github/pato.gif"  height="64pt" style="margin-bottom: -28%"> PATO: high PerformAnce TriplexatOr

![gnu workflow](https://img.shields.io/github/actions/workflow/status/UDC-GAC/pato/compile-and-test-gnu.yml?label=gnu)
![clang workflow](https://img.shields.io/github/actions/workflow/status/UDC-GAC/pato/compile-and-test-clang.yml?label=clang)
[![issues](https://img.shields.io/github/issues/UDC-GAC/pato)](https://github.com/UDC-GAC/pato/issues)
[![doi](https://img.shields.io/badge/doi-btad134-blue)](https://doi.org/10.1093/bioinformatics/btad134)
[![release](https://img.shields.io/github/release/UDC-GAC/pato?include_prereleases=&sort=semver&color=blue)](https://github.com/UDC-GAC/pato/releases/)
[![license](https://img.shields.io/badge/license-MIT-blue)](LICENSE)
![platform](https://img.shields.io/badge/platform-linux%20|%20linux--arm%20|%20macos%20|%20macos--m1%20|%20windows%20(wsl)-blue)

PATO: high PerformAnce TriplexatOr is a high performance tool for the fast and efficient detection of triple helices and triplex features in nucleotide sequences. PATO: high PerformAnce TriplexatOr is a modern alternative to [Triplexator](https://github.com/Gurado/triplexator) and [TDF](https://github.com/CostaLab/reg-gen) and functions as a drop in replacement to accelerate the triplex analyses in multicore computers. Moreover, PATO: high PerformAnce TriplexatOr's efficiency allows a more exhaustive search of the triplex-forming solution space, and so it achieves higher levels of prediction accuracy in far less time than other tools in the state of the art.

## Version
Version v0.0.2.

PATO: high PerformAnce TriplexatOr v0.0.2 has been tested and is fully functional on the following platforms:
* Linux x86_64
* Linux ARMv8
* Windows x86_64 (using [WSL](https://learn.microsoft.com/en-us/windows/wsl/install))
* macOS x86_64
* macOS M1

## Requirements
PATO: high PerformAnce TriplexatOr does not have any hardware requirements, and it can run on, but is not limited to, x86, x86_64, and ARM processors.

To compile and execute PATO: high PerformAnce TriplexatOr, the following software is required:
* GNU Make
* C++ compiler ([c++17 compliant](https://en.cppreference.com/w/cpp/compiler_support/17) and with support for OpenMP directives)

A valid combination of these tools may be: GNU Make v3.82, and GCC v9.3.0.

### macOS requirements
Because Apple's `clang` compiler does not support OpenMP parallelism, macOS users are required to install an OpenMP-enabled compiler in order to build PATO: high PerformAnce TriplexatOr from source. To do so, we recommend installing the GNU Compiler Collection v12 (or any other version with support for the c++17 standard) using Homebrew as follows:
```bash
brew install gcc@12
```

Afterwards, make sure `gcc-12` and `g++-12` are included in your system's `PATH` variable by executing:
```bash
echo "export PATH=$(brew --prefix gcc@12)/bin:\$PATH" >> ~/.zshrc && source ~/.zshrc
```

Lastly, check that you now indeed have a working Homebrew `gcc-12` and `g++-12` installation by executing:
```bash
gcc-12 --version
gcc-12 (Homebrew GCC 12.2.0) 12.2.0
Copyright (C) 2022 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

and

```bash
g++-12 --version
g++-12 (Homebrew GCC 12.2.0) 12.2.0
Copyright (C) 2022 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

## Compilation
Download the source code from this repository, either use Git or download a copy from GitHub, and let GNU Make automatically compile PATO: high PerformAnce TriplexatOr for you:
```bash
make gnu -j$(getconf _NPROCESSORS_ONLN)
```

### macOS compilation
macOS users must explicitly specify an OpenMP-enabled compiler to compile PATO: high PerformAnce TriplexatOr. For instance, to use the Homebrew `g++-12` compiler, execute:
```bash
make gnu -j$(getconf _NPROCESSORS_ONLN) CXX=g++-12 LD=g++-12
```

## Execution
Now that PATO: high PerformAnce TriplexatOr has been compiled, execute the application as follows:
```bash
./target/gnu/pato.release [options] {-ss tfo_file | -ds tts_file | -ss tfo_file -ds tts_file}
```

Execute `./target/gnu/pato.release --help` for a detailed list of execution modes, options, and flags.

### Setting the number of threads
PATO: high PerformAnce TriplexatOr uses OpenMP to parallelize its triplex search algorithms. The OpenMP runtime will automatically spawn as many threads as there are available CPU cores. To reduce the number of threads spawned by the application, therefore, one has to explicitly set the `OMP_NUM_THREADS` environment variable to a value greater than 0. For instance, to run PATO: high PerformAnce TriplexatOr with 4 threads, execute:
```bash
OMP_NUM_THREADS=4 ./target/gnu/pato.release ...
```

### Setting the number of simultaneous sequences
To reduce the memory footprint of PATO: high PerformAnce TriplexatOr, one can set the maximum number of sequences that may processed simultaneously by the triplex search algorithms. This is done by setting the `-cs` or `--chunk-size` option to a value greater than 0 (128 by default). For instance, to process a dataset in chunks of 32 sequences, execute:
```bash
./target/gnu/pato.release --chunk-size 32 ...
```

To give an upper bound of the memory consumption of PATO: high PerformAnce TriplexatOr, one can use the following formula:

$$
\begin{flalign}
& \text{mem}(cs, t, l) = cs \cdot \text{tpx}(l) + cs \cdot \text{len}(l) + \begin{cases}
t \cdot \text{len}(l), & \text{if } cs > t\\
cs \cdot \text{len}(l), & \text{if } cs \leq t
\end{cases} &
\end{flalign}
$$

where $cs$ is the number of simultaneous sequences, $t$ is the number of threads, $l$ is the longest sequence in a given dataset, and $\text{tpx}(l)$ is the size of the triplex features found in $l$.

It is possible to further reduce the memory usage of the application by disabling the filtering of low-complex regions. This can be done by setting the `-fr` or `--filter-repeats` option to `false`. In such a case, sequences should be filtered before being passed to PATO: high PerformAnce TriplexatOr (the [Ensembl genome browser](https://www.ensembl.org) provides filtered sequences) and the formula becomes:

$$
\begin{flalign}
& \text{mem}(cs, t, l) = cs \cdot \text{tpx}(l) + cs \cdot \text{len}(l) &
\end{flalign}
$$

In general, one can't go wrong by setting the number of simultaneous sequences to a value equal to the number of threads that PATO: high PerformAnce TriplexatOr is going to use. However, if the sequences of a dataset are very long, it may be necessary to reduce the number of simultaneous sequences to avoid running out of memory.

If you are unsure about the number of simultaneous sequences to use, you can set the `-cs` or `--chunk-size` option to 1. Although this may hurt parallelism by a small amount, it will allow you to run PATO: high PerformAnce TriplexatOr on any dataset without having to worry about the memory footprint of the application.

## Test
If you are interested in testing the correctness of your installation, run the following commands:
```bash
make gnu -j$(getconf _NPROCESSORS_ONLN) BUILD=release
make gnu -j$(getconf _NPROCESSORS_ONLN) BUILD=serial
bash test/test.bash gnu
```

## Cite us
If you use PATO: high PerformAnce TriplexatOr in your research, please cite our work using the following reference:
```
Iñaki Amatria-Barral, Jorge González-Domínguez and Juan Touriño. "PATO: Genome-wide prediction of lncRNA-DNA triple helices". Bioinformatics (2023)
```

## License
PATO: high PerformAnce TriplexatOr is free software and as such it is distributed under the [MIT License](LICENSE). However, PATO: high PerformAnce TriplexatOr makes use of several modules which are not original pieces of work. Therefore, their usage is subject to their corresponding [THIRDPARTYLICENSE](THIRDPARTYLICENSES) and all rights are reserved to their authors.
