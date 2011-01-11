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

use warnings;

use strict;

die "Arguments invalid" if @ARGV != 2;

my $ped= shift @ARGV;
my $gen= shift @ARGV;

open (PED, "<$ped") or die "Impossible to read from the pedigree file $ped.";
open (GEN, "<$gen") or die "Impossible to read from the genotype file $gen.";

$_= <GEN>;
$_ =~ /([0-9]+) ([0-9]+)/ or die "Impossible to interpret the header $_.";

my $n_gen= $1;
my $n_loci= $2;

print "# no. loci\n";
print "$n_loci\n";
print "# no. individuals\n";
print "$n_gen\n";

print "# genotypes\n";

while (<GEN>) {
    s/ //g;
    print "$_";
}

print "# pedigree\n";

$_= <PED>;
while (<PED>) {
    chomp;
    next if /N\tN/;
    next if /$^/;
    /([0-9]+)\t[MF]\t([0-9]+)\t([0-9]+)/ or next;
    my $id= $1-1;
    my $fid= $2-1;
    my $mid= $3-1;
    print "$fid $mid $id\n";
}
print "# END\n";
