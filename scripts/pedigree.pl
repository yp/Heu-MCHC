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

my %genotypes=();
my $count= 0;

while (<>) {
    last if /^# ?genotypes/;
}

print STDERR "digraph G{\n";
print STDERR "  node [shape=record];\n";
print STDERR "  graph [ordering=out,rankdir=BT,nodesep=2];\n";
print STDERR "  edge [minlen=2,dir=back];\n";
#Read the genotypes
while (<>) {
    chmod;
    last if m/^# ?pedigree/i;
    my $g= $_;
    my @gen =();
    print STDERR "  n$count [label=\"{G_$count";
    my $i= 0;
    while ($g =~ s/^(.)//) {
        print STDERR "|";
        push(@gen, $1);
        print STDERR "00" if ($1 eq "0");
        print STDERR "11" if ($1 eq "1");
        print STDERR "01" if ($1 eq "2");
        $i++;
    }
    print STDERR "}\"];\n";
    $genotypes{$count}= \@gen;
    $count++;
}
print "Genotypes read: $count\n";

#Read the trios
while (<>) {
    print $_;
    last if m/^# END/i;
    m/([0-9]+) ([0-9]+) ([0-9]+)/ or next;
    print STDERR "  n$3->n$1 [label=\"F\"];\n";
    print STDERR "  n$3->n$2 [label=\"M\"];\n";
}
print STDERR "}\n";
