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
## summary-pedigree.sh
## 
## Made by Yuri
## Login   <yuri@dottxxii-12.dottorato.disco.unimib.it>
## 
## Started on  Fri Aug 28 10:57:35 2009 Yuri
## Last update Fri Aug 28 10:57:35 2009 Yuri
##

f=`basename $1`;
cart=`dirname $1`;
suff=$3

header=${2:-0};

if test $header -eq 0; then
    echo -ne "dir\t ped\t indiv\t loci\t mrate\t rrate\t addtl_child_prob\t multimat_prob\t";
    echo -ne "sim_mut\t sim_rec\t sim_ev\t min_mut\t min_rec\t min_ev\t pred_mut\t pred_rec\t pred_ev\t";
    echo -e "common_ev\t only_sim_ev\t only_pred_ev\t tot_loci\t tot_het_loci\t err_phases\t time_centisec"
fi

echo -ne "${cart}/\t${f}\t"

info-from-pedigree.sh ${cart}/$f;

for ev in ${cart}/events/min-out-$f ${cart}/pred/out-${suff}$f; do
    gawk '/= MUT/ { ++mut } /= REC/ { ++rec } END { printf "\t%d\t%d\t%d", mut, rec, mut+rec }' $ev;
done

echo -ne "\t";

compare-events.pl ${cart}/events/min-out-$f ${cart}/pred/out-${suff}$f 2> /dev/null;

compare-results.pl ${cart}/events/min-hap-out-$f ${cart}/pred/hap-${suff}out-$f | gawk '
{ suml += $2; sum2 += $3; sumd += $4; if ($5 != $4) print "****ERROR****"; }
END { printf "\t%d\t%d\t%d", suml, sum2, sumd }'

echo -ne "\t"

tail ${cart}/pred/log-${suff}$f | grep "timer general" | sed 's/ \+/ /g' | cut -d ' ' -f 6 | sed 's/[0-9][0-9][0-9][0-9]$//'

