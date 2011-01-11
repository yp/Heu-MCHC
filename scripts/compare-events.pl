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

@ARGV >= 2 or die "Please specify two parameters: FILE-EVENTS-1 FILE-EVENTS-2\n";

my $file1= $ARGV[0];
my $file2= $ARGV[1];

open F1, "<$file1" or die "I cannot open file $file1\n";
open F2, "<$file2" or die "I cannot open file $file2\n";

my %ev1= ();
my %ev2= ();

while (<F1>) {
    chomp;
    m/^= (MUT|REC)/ or next;
    s/^= //;
    s/[^A-Za-z0-9]+/_/g;
    $ev1{$_}= $_;
}

while (<F2>) {
    chomp;
    m/^= (MUT|REC)/ or next;
    s/^= //;
    s/[^A-Za-z0-9]+/_/g;
    $ev2{$_}= $_;
}

my $common= 0;
my $only1= 0;
my $only2= 0;
my $sonly1= "";
foreach my $e1 (keys %ev1) {
    if (exists $ev2{$e1}) {
        print STDERR "common $e1\n";
        ++$common;
    } else {
        $sonly1 .= "ev1 $e1\n";
        ++$only1;
    }
}
print STDERR "$sonly1";
foreach my $e2 (keys %ev2) {
    if (exists $ev1{$e2}) {
    } else {
        print STDERR "ev2 $e2\n";
        ++$only2;
    }
}
print "$common\t$only1\t$only2";

