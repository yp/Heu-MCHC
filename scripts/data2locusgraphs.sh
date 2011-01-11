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

#echo "digraph G {" > tot-pedigree.dot;
for l in `seq 1 15`; do
   cat pedigree.dot | grep -v "> f" | sed 's/f\(...\)m\(...\) -> i\(...\)/i\1 -> i\3;i\2 -> i\3/' | grep -v "f" | sed 's/i0\+\([0-9]\)/i\1/g' > tmp.dot;
   grep "G " gt-ped-conf-426.txt | sed 's/.\+G //' | cut -b $l | gawk '{print "i" NR-1, $0}' > tmp.locus;
   while read i g; do
	if test $g -ne 2; then
	sed "s/$i\[/$i\[label=\"$i($g)\",/" tmp.dot  | sed "s/$i\([ ;\[]\)/${i}L$l\1/g" > tmp2.dot;
	else
	sed "s/$i\[/$i\[label=\"$i($g)\",color=gray,style=\"filled\"/" tmp.dot  | sed "s/$i\([ ;\[]\)/${i}L$l\1/g" > tmp2.dot;
	fi;
	mv tmp2.dot tmp.dot;
   done < tmp.locus;
   cat tmp.dot | sed "s/digraph G/digraph cluster$((l-1))/" | sed "s/}/label=\"locus $((l-1))\"; }/" >> tot-pedigree.dot;
done
#echo "}" >> tot-pedigree.dot
rm tmp.dot tmp.locus




