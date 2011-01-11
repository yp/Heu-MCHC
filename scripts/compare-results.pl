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

sub hd {
    return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}

sub mydiff {
    my ($rt1, $rt2)= @_;
    my %diffs= ();
    if (join(":", sort keys %$rt1) eq join(":", sort keys %$rt2)) {
        foreach my $id (sort {$a <=> $b} keys %$rt1) {
            $diffs{$id}= hd($rt1->{$id},$rt2->{$id});
        }
    } else {
        die "The lists of individuals are different.\n".
            "Read L1=" . join(":", sort keys %$rt1) . "\n" .
            "Read L2=" . join(":", sort keys %$rt2) . "\n";
    }
    return \%diffs;
}

@ARGV >= 2 or die "Please specify two parameters: FILE-RIS-1 FILE-RIS-2\n";

my $file1= $ARGV[0];
my $file2= $ARGV[1];

my $xorout= 1;
if (@ARGV>=3) {
    $xorout= $ARGV[2];
}

open F1, "<$file1" or die "I cannot open file $file1 (the other file is $file2)\n";
open F2, "<$file2" or die "I cannot open file $file2 (the other file is $file1)\n";

my %g1= ();
my %p1= ();
my %m1= ();
my %g2= ();
my %p2= ();
my %m2= ();
my %f= ();

while (<F1>) {
    chomp;
    s/^ +//;
    my @a= split " ";
    if ($a[1] eq "FOUNDER") {
        $f{$a[0]}= 1;
    } elsif ($a[1] eq "G") {
        $g1{$a[0]}= $a[2];
    } elsif ($a[1] eq "P") {
        $p1{$a[0]}= $a[2];
    } elsif ($a[1] eq "F") {
        if (defined $m1{$a[0]}) {
            $p1{$a[0]}= $m1{$a[0]};
        } 
        $m1{$a[0]}= $a[2];
    } elsif ($a[1] eq "M") {
        if (exists $m1{$a[0]}) {
            $p1{$a[0]}= $a[2];
        } else {
            $m1{$a[0]}= $a[2];
        }
    }
}
while (<F2>) {
    chomp;
    s/^ +//;
    my @a= split " ";
    if ($a[1] eq "FOUNDER") {
        $f{$a[0]}= 1;
    } elsif ($a[1] eq "G") {
        $g2{$a[0]}= $a[2];
    } elsif ($a[1] eq "P") {
        $p2{$a[0]}= $a[2];
    } elsif ($a[1] eq "F") {
        if (exists $m2{$a[0]}) {
            $p2{$a[0]}= $m2{$a[0]};
        } 
        $m2{$a[0]}= $a[2];
    } elsif ($a[1] eq "M") {
        if (exists $m2{$a[0]}) {
            $p2{$a[0]}= $a[2];
        } else {
            $m2{$a[0]}= $a[2];
        }
    }
}

# Check genotypes differences
my $gd= 0;
my $gdr= mydiff(\%g1, \%g2, $file1, $file2);
$gd += $_ foreach values %$gdr;

$gd == 0 or die "The genotype set contains differences\n. Differences ".join(" ",values %$gdr)."\n";

# First heterozygous locus of a founder is 0 in paternal haplotype
 for my $id (sort {$a <=> $b} keys %f) {
     if (hd($p1{$id}, $p2{$id})>hd($p1{$id}, $m2{$id})) {
         my $t= $p1{$id};
         $p1{$id}= $m1{$id};
         $m1{$id}= $t;
     }
 }

my $pd= 0;
my $pdr= mydiff(\%p1, \%p2, $file1, $file2);
$pd += $_ foreach values %$pdr;


my $md= 0;
my $mdr= mydiff(\%m1, \%m2, $file1, $file2);
$md += $_ foreach values %$mdr;

(join(":", sort keys %$pdr) eq join(":", sort keys %$mdr)) or die "The haplotype sets refer to different individuals.\n";

foreach my $id (sort { $a <=> $b } keys %$pdr) {
    print "$id\t" . $g1{$id} =~ tr/012// . "\t" . $g1{$id} =~ tr/2// .
        "\t$pdr->{$id}\t$mdr->{$id}\t". ($pdr->{$id}+$mdr->{$id})."\t";
    my $xor= $p1{$id} ^ $p2{$id};
    $xor =~ tr/\000\001-\255/.*/;
    print $p1{$id} . " " . $p2{$id} . " $xor\t";
    $xor= $m1{$id} ^ $m2{$id};
    $xor =~ tr/\000\001-\255/.*/;
    print $m1{$id} . " " . $m2{$id} . " $xor" if $xorout==1;
    print "\n";
}

