
   Heu-MCHC
==============

A fast and accurate heuristic algorithm for the haplotype inference
problem on pedigree data with recombinations and mutations

by  [Yuri Pirola](http://bimib.disco.unimib.it/index.php/Pirola_Yuri)
and [Tao Jiang](http://www.cs.ucr.edu/~jiang/)

------------------------------------------------------------------------

## Introduction ##

Heu-MCHC is a fast and accurate heuristic for the Minimum-Change
Haplotype Configuration (MCHC) problem, i.e. a combinatorial formulation
of the haplotype inference problem on pedigree data where the total
number of recombinations and point mutations has to be minimized.


## Reference ##

Please cite the following paper when using Heu-MCHC.

Yuri Pirola, Paola Bonizzoni, and Tao Jiang.
Haplotype Inference on Pedigrees with Recombinations and Mutations.
In Algorithms in Bioinformatics, 10th Int. Workshop, WABI 2010, Liverpool UK, Sep 6-8, 2010, Proceedings. V. Moulton and M. Singh (Ed.).
Springer. Vol. 6293 pp. 148-161 (2010).
[Link](http://dx.doi.org/10.1007/978-3-642-15294-8_13)


## License ##

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.


## General Overview ##

Heu-MCHC is a software pipeline composed by two distinct steps:

*  ped-hi, which, given a genotyped pedigree, predicts a set of
   variation events (recombinations and point mutations) that allows to
   recover a valid haplotype configuration.
*  to-hap, which, given a genotyped pedigree and the set of variation
   events, generates a haplotype configuration that is consistent with
   the genotyped pedigree and that induces the set of variation events
   given as input.


### Usage ###

To infer haplotypes starting from a genotyped pedigrees, both programs
have to be used.
First, we must invoke `ped-hi` which reads from the standard input the
genotyped pedigree (please refer below for the description of the file
format) and writes the predicted variation events to a file whose name
has been specified as first command-line parameter.
During the execution, `ped-hi` prints to standard error a description of
the steps it performs.
For example, with the following command line

    ./ped-hi out-events.txt < genotyped-pedigree.txt 2> log-ped-hi.txt

the program reads the genotyped pedigree from file
`genotyped-pedigree.txt`, writes the predicted events to the file
`out-events.txt`, and logs its activities to the file `log-ped-hi.txt`.

Second, we must invoke `to-hap` in order to recover a set of haplotypes
which induces the predicted variation events.
Program `to-hap` requires two parameters specified by means of the two
flags `-p` and `-e`.
Flag `-p` specifies the name of the file which describes the genotyped
pedigree, while flag `-e` specifies the name of the file which contains
the set of predicted variation events. To read one of the two files from
standard input, we must specify `-` as filename.
A consistent haplotype configuration is printed to standard output if
the input data are valid, otherwise it terminates with an error.
Like `ped-hi`, `to-hap` prints to standard error a description of the
steps it performs.

For example, if the genotyped pedigree is described in file
`genotyped-pedigree.txt`, the following commands perform the complete
haplotype inference process.

    ./ped-hi out-events.txt < genotyped-pedigree.txt 2> log-ped-hi.txt
    ./to-hap -p genotyped-pedigree.txt -e out-events.txt > out-haplotypes.txt 2> log-to-hap.txt


### File formats ###

#### Genotyped Pedigrees ####
The genotyped pedigree is described by a single file, which contains
both the pedigree structure and the set of genotypes.
We describe the file format by a simple example.

    2
    7
    #genotypes
    10
    02
    11
    22
    22
    22
    22
    #pedigree
    0 1 4
    2 3 5
    4 5 6
    #end

The first row specifies the number of loci (2, in this case) over which
the genotypes are defined.
The second row (7, in this case) specifies the pedigree size (i.e., the
number of its members).
Rows starting with `#` are considered as comments and are discarded.
Genotypes are specified from the third row (discarding comments).
Each row represents the multi-locus genotype of a single individual and
genotypes are encoded as follows: 0=major allele homozygous, 1=minor
allele homozygous, 2=heterozygous.
The first genotype (the one in third row) is associated to individual
0, the second genotype (fourth row) is associated to individual 1, and
so on.
After all the genotypes have been described, the file specifies the
pedigree structure by means of a series of triplets (each one on a
single row).
Each triplet specifies a single trio, where the first number is the
identifier of the father, the second one is the identifier of the
mother, and the third one is the identifier of the child.
Identifiers start from 0 and are associated to genotypes as described
above.

#### Set of Variation Events ####
A set of variation events is described by a file containing a series of
rows with the following structure:

    = VAR i j k

where:

- =   is a fixed character
- VAR is an abbreviation of the type of variation event (REC for
   recombinations and MUT for point mutations),
- i   is the identifier of the individual where the event has
   finally occurred,
- j   is the identifier of the locus where the event has occurred
   (0-based),
- k   is a flag which specifies if the event has occurred between
   individual i and his father (k=0) or between i and his mother (k=1).

A row which does not conform to such a structure is silently ignored.

### Haplotype configuration ###
A haplotype configuration is described by a space-separated file in
which rows have the following structure:

    i T "remaining"

where:

- i  is the identifier of the individual which the row refers to
- T  is a flag which specifies the content of the remaining part of
   the row. If T='G', then "remaining" is the genotype f the
   individual i; if T='P', then "remaining" is the paternal
   haplotype of the individual; if T='M', then "remaining" is the
   maternal haplotype of the individual; if T='FOUNDER', then the
   individual is a pedigree founder (hence its haplotypes might
   be swapped, since there is no way to distinguish the paternal
   and the maternal haplotype in a pedigree founder).
-  "remaining"  is the remaining part of the row, and its meaning is
   determined as described above.


## Package Contents ##

A binary distribution is currently provided for your convenience at the
[AlgoLab web site](http://www.algolab.eu/Heu-MCHC/).
The software is available only for Linux 64-bit operating systems (fully
tested on Ubuntu 10.04 x64).

The distribution provides an implementation of the two steps in
different `flavors`, denoted by a suffix of the program name.
Below we list the flavors currently available:

- `regular`/`fast`: in the `fast` flavor the Sum-Product algorithm is
  iterated 3 times instead of the 10 times of the `regular` flavor. As
  the name suggest, the `fast` flavor is faster than the `regular` one
  but it might be also less accurate.
- `complete`/`only_rec`/`only_mut`: the `complete` flavor allows both
  recombinations and mutations, the `only_rec` flavor allows only
  recombination events, while the `only_mut` flavor allows only point
  mutations.
- `base`/`fw`: the `base` flavor implements the Sum-Product algorithm
  from scratch, while the `fw` flavor uses the Framewave mathematical
  library. Please notice that, in order to use the `fw` flavor, the
  Framewave library must be installed in the system. The library can be
  freely downloaded from [SourceForge](http://framewave.sourceforge.net/)
  and it is also available from the Debian/Ubuntu repositories.
  The `fw` flavor is functionally identical to the `base` flavor but,
  in general, is considerably faster (especially for larger instances).
  We strongly suggest to use the `fw` flavor, whenever it is possible.
  Notice that different results might be obtained in exceptional
  circumstances due to different rounding and accuracy of intermediate
  results.

We also provide various simple genotyped pedigrees as example (in the
`examples` directory).


## Contacts ##

Please contact Yuri Pirola for additional information.  
E-mail:   `yuri.pirola (*) gmail.com`  
Web page: `http://bimib.disco.unimib.it/index.php/Pirola_Yuri`
