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
## exec-pedphase.sh
## 
## Made by Yuri
## Login   <yuri@dottxxii-12.dottorato.disco.unimib.it>
## 
## Started on  Wed Jul 29 15:57:57 2009 Yuri
## Last update Wed Jul 29 15:57:57 2009 Yuri
##


## Given the same input as ped-hi, convert the data, execute PedPhase and postprocess the output

fin=$1

dir=$2

if test ! -r $fin; then
    echo "Impossible to find input file.";
    exit -1;
fi;

${dir}ped2pedphase < $fin > pp-input-$fin
{
time wine ${dir}PedPhase.exe -i pp-input-$fin > log-pp-$fin 2>&1
} 2> log-time-pp-$fin 

mv pp-input-`echo $fin | sed 's/...$/out/'` rawout-pp-$fin
rm pp-input-`echo $fin | sed 's/...$/dat/'`
rm pp-input-`echo $fin | sed 's/...$/res/'`
rm pp-input-`echo $fin | sed 's/...$/msp/'`
rm PedPhaseILP.log

${dir}pedphase2hap.pl < rawout-pp-$fin > hap-pp-$fin
${dir}pedphase2rec.gawk < rawout-pp-$fin > out-pp-$fin




