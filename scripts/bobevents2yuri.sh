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
##
## bobevents2yuri.sh
## 
## Made by Yuri
## Login   <yuri@dottxxii-12.dottorato.disco.unimib.it>
## 
## Started on  Sat Aug 29 00:48:28 2009 Yuri
## Last update Sat Aug 29 00:48:28 2009 Yuri
##


for p in $*; do 
f=`basename $p`;
d=`dirname $p`;
echo "= SIMULATED EVENTS" > $d/events/events-$f;
gawk 'BEGIN {FS="_"}
/Recombinants:/ {p=0} 
// {if (p==1 && NF>0) {par=$1-1; chi=$2-1; locus=$3+0; printf "%d\t%d\t%d\n", par, chi, locus;}} 
/Mutations:/ {p=1}' $d/events/${f/ped-gen/hap} | 
  while read par chi loc; do
  printf "= MUT %4d %4d  " ${chi} ${loc} >> $d/events/events-$f;
  tail -n+$((`grep -n "# pedigree" $p | sed 's/:.\+//'`+1)) $p | 
    head -n-1 | 
    grep "^[0-9]\+ $par $chi" > /dev/null && echo "1" >> $d/events/events-$f || echo "0" >> $d/events/events-$f;
  done;
done
