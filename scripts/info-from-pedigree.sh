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
# along with Heu-MCHC.  If not, see <http://www.gnu.org/licenses/>.
#
####
##
## info-from-pedigree.sh
## 
## Made by Yuri
## Login   <yuri@dottxxii-12.dottorato.disco.unimib.it>
## 
## Started on  Mon Jul 27 11:40:21 2009 Yuri
## Last update Mon Jul 27 11:40:21 2009 Yuri
##

f=${1:-pedigree.txt}

if test ! -f $f; then
    echo "Pedigree not found!";
    exit;
fi;

head $f | sed 's/ [ ]\+/ /g' | cut -d ' ' -f 3 | {
read h; # header UNUSED
read conf; # configuration number UNUSED
read pmut;
read prec;
read nind;
read nloc;
read pchildren;
read pmulti;
echo -ne "$nind\t$nloc\t$pmut\t$prec\t$pchildren\t$pmulti\t";
}

echo -ne "`grep "^#\* MUT" $f | wc -l`\t";
echo -ne "`grep "^#\* REC" $f | wc -l`\t";
echo -ne "`grep "^#\*" $f | wc -l`";


