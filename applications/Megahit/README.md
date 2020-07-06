MEGAHIT Benchmarks
=======

Installation
---------------

### Building from source

#### Prerequisites

-   For building: zlib, cmake &gt;= 2.8, g++ &gt;= 4.8.4
-   For running: gzip and bzip2

```sh
./install.sh
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release  # add -DCMAKE_INSTALL_PREFIX=MY_PREFIX if needed
make benchmarks -j4
make test_load_reads
make test_load_seq
make test_build_sdbg
make test_assemble
```
