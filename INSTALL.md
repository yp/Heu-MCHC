
   Heu-MCHC
==============

A fast and accurate heuristic algorithm for the haplotype inference
problem on pedigree data with recombinations and mutations

------------------------------------------------------------------------


## Compilation ##

The compilation process is divided in two steps:

1. Compilation and installation of the `m4ri` library
2. Compilation of the program `Heu-MCHC`


### Compilation and installation of the `m4ri` library ###

In this project, we use a slightly modified version of the `m4ri`
library which is included in the `pers-lib/m4ri/` directory.
Assuming that the required dependencies are satisfied, it suffices to
use the following commands:

    cd pers-lib/m4ri
    ./configure --prefix=${PWD}/../local/
    make && make install
    cd ../..

### Compilation of the program `Heu-MCHC` ###

As explained in the `README` file, `Heu-MCHC` can be compiled in
different flavors by setting some environmental variables during the
invocation of the `make` process.
The most common combination of flavors can be produced by issuing
the following command:

    make STATUS=production LINKING=static EXT_DEFINE="-DMAX_TENT=10" MATH_LIB=none RECOMBINATIONS=yes MUTATIONS=yes reall

A faster alternative that uses the `FRAMEWAVE` library can be produced
by the following command:

    make STATUS=production LINKING=static EXT_DEFINE="-DMAX_TENT=10" MATH_LIB=framewave RECOMBINATIONS=yes MUTATIONS=yes reall


The various environmental variables and their values are quite
self-explicative: please refer to the `Makefile` for further
information.




## Automatic Compilation of All Flavors ##

The simple script `./scripts/deploy.sh` automatically creates a
directory containing the binaries compiled for every interesting
combination of flavors (as described in the `README` file).


