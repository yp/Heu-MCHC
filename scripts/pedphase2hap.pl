#!/usr/bin/perl -w
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

use strict;

my $previd= -1;
while (<>) {
    chomp;
    my @arr= split /\t/;
    $#arr -= 2;
    my $id= $arr[1]-1;
    last if $id < $previd;
    $previd= $id;
    my $pid= $arr[2];
    my $mid= $arr[3];
    if ($pid == 0 or $mid == 0 ) {
        ($mid+$pid) == 0 or die "Individual $id is a partial founder.\n";
        printf "%5d FOUNDER\n", $id;
    }
    my $g= "";
    my $p= "";
    my $m= "";
    for (my $i= 7; $i<=$#arr; ++$i) {
        $arr[$i] =~ /^([1-2])\|([1-2])$/ or die "Allele not recognized!\n";
        if ($1 eq $2) {
            $g = $g . ($1-1);
        } else {
            $g= $g. "2";
        }
        $p .= ($1-1);
        $m .= ($2-1);
    }
    printf "%5d G $g\n", $id;
    printf "%5d P $p\n", $id;
    printf "%5d M $m\n", $id;
}
