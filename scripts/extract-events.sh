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
## extract-events.sh
## 
## Made by Yuri
## Login   <yuri@dottxxii-12.dottorato.disco.unimib.it>
## 
## Started on  Thu Aug 27 14:58:56 2009 Yuri
## Last update Thu Aug 27 14:58:56 2009 Yuri
##

mkdir -p events
for f in $*; do
    echo -ne "$f\t";
    intest=`grep "^#\*" $f | gawk '/#\* MUT / {mut++} /#\* REC / {rec++} END { printf "mut %3d rec %3d  tot %3d\n", mut, rec, mut+rec;}'`;
    echo "= SIMULATED EVENTS  $intest" > events/events-$f;
    grep "^#\*" $f | sed 's/#\*/=/' >> events/events-$f;
    echo "$intest";
done

