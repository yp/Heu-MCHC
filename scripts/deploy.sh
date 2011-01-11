#!/bin/bash
####
#
#
#                               Heu-MCHC
#
# A fast and accurate heuristic algorithm for the haplotype inference
# problem on pedigree data with recombinations and mutations
#
# Copyright (C) 2009,2010,2011  Yuri PIROLA
#
# Distributed under the terms of the GNU General Public License (GPL)
#
#
# This file is part of Heu-MCHC.
#
# Heu-MCHC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Heu-MCHC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
#
####
LANG=C

DEST=Heu-MCHC

rm -rf ${DEST}

mkdir ${DEST}

cp -r docs/* ${DEST}/

echo "Package:       Heu-MCHC" > ${DEST}/version.txt
echo "Code revision: `git branch -v | grep '^\*' | awk '{print $3}'`" >> ${DEST}/version.txt
echo "Build date:    `date`" >> ${DEST}/version.txt
echo "Author:        `git config --get user.name`" >> ${DEST}/version.txt


make STATUS=production LINKING=static BITNESS=32 EXT_DEFINE="-DMAX_TENT=3" MATH_LIB=none RECOMBINATIONS=yes MUTATIONS=yes reall
DESC=complete-fast-base
cp bin/ped-hi ${DEST}/ped-hi-${DESC}
cp bin/to-hap ${DEST}/to-hap-${DESC}

make STATUS=production LINKING=static BITNESS=32 EXT_DEFINE="-DMAX_TENT=3" MATH_LIB=none RECOMBINATIONS=yes MUTATIONS=no reall
DESC=only_rec-fast-base
cp bin/ped-hi ${DEST}/ped-hi-${DESC}
cp bin/to-hap ${DEST}/to-hap-${DESC}

make STATUS=production LINKING=static BITNESS=32 EXT_DEFINE="-DMAX_TENT=3" MATH_LIB=none RECOMBINATIONS=no MUTATIONS=yes reall
DESC=only_mut-fast-base
cp bin/ped-hi ${DEST}/ped-hi-${DESC}
cp bin/to-hap ${DEST}/to-hap-${DESC}

make STATUS=production LINKING=static BITNESS=32 EXT_DEFINE="-DMAX_TENT=3" MATH_LIB=framewave RECOMBINATIONS=yes MUTATIONS=yes reall
DESC=complete-fast-fw
cp bin/ped-hi ${DEST}/ped-hi-${DESC}
cp bin/to-hap ${DEST}/to-hap-${DESC}

make STATUS=production LINKING=static BITNESS=32 EXT_DEFINE="-DMAX_TENT=3" MATH_LIB=framewave RECOMBINATIONS=yes MUTATIONS=no reall
DESC=only_rec-fast-fw
cp bin/ped-hi ${DEST}/ped-hi-${DESC}
cp bin/to-hap ${DEST}/to-hap-${DESC}

make STATUS=production LINKING=static BITNESS=32 EXT_DEFINE="-DMAX_TENT=3" MATH_LIB=framewave RECOMBINATIONS=no MUTATIONS=yes reall
DESC=only_mut-fast-fw
cp bin/ped-hi ${DEST}/ped-hi-${DESC}
cp bin/to-hap ${DEST}/to-hap-${DESC}


make STATUS=production LINKING=static BITNESS=32 EXT_DEFINE="-DMAX_TENT=10" MATH_LIB=none RECOMBINATIONS=yes MUTATIONS=yes reall
DESC=complete-regular-base
cp bin/ped-hi ${DEST}/ped-hi-${DESC}
cp bin/to-hap ${DEST}/to-hap-${DESC}

make STATUS=production LINKING=static BITNESS=32 EXT_DEFINE="-DMAX_TENT=10" MATH_LIB=none RECOMBINATIONS=yes MUTATIONS=no reall
DESC=only_rec-regular-base
cp bin/ped-hi ${DEST}/ped-hi-${DESC}
cp bin/to-hap ${DEST}/to-hap-${DESC}

make STATUS=production LINKING=static BITNESS=32 EXT_DEFINE="-DMAX_TENT=10" MATH_LIB=none RECOMBINATIONS=no MUTATIONS=yes reall
DESC=only_mut-regular-base
cp bin/ped-hi ${DEST}/ped-hi-${DESC}
cp bin/to-hap ${DEST}/to-hap-${DESC}

make STATUS=production LINKING=static BITNESS=32 EXT_DEFINE="-DMAX_TENT=10" MATH_LIB=framewave RECOMBINATIONS=yes MUTATIONS=yes reall
DESC=complete-regular-fw
cp bin/ped-hi ${DEST}/ped-hi-${DESC}
cp bin/to-hap ${DEST}/to-hap-${DESC}

make STATUS=production LINKING=static BITNESS=32 EXT_DEFINE="-DMAX_TENT=10" MATH_LIB=framewave RECOMBINATIONS=yes MUTATIONS=no reall
DESC=only_rec-regular-fw
cp bin/ped-hi ${DEST}/ped-hi-${DESC}
cp bin/to-hap ${DEST}/to-hap-${DESC}

make STATUS=production LINKING=static BITNESS=32 EXT_DEFINE="-DMAX_TENT=10" MATH_LIB=framewave RECOMBINATIONS=no MUTATIONS=yes reall
DESC=only_mut-regular-fw
cp bin/ped-hi ${DEST}/ped-hi-${DESC}
cp bin/to-hap ${DEST}/to-hap-${DESC}

pushd ${DEST}
ln -s ped-hi-complete-regular-base ped-hi
ln -s to-hap-complete-regular-base to-hap
popd

tar cvjf ${DEST}.tar.bz2 ${DEST}/
tar cvzf ${DEST}.tar.gz ${DEST}/

