#!/bin/perl

open(F, $ARGV[0]) || die "Can't open file\n";

while($F = <F>) {

  chomp($F);

  ($q, $allele) = split(/\t/, $F);

  $Hash{$q}{$allele} ++

}

foreach $gene (keys %Hash) {

  foreach $allele (keys %{$Hash{$gene}}) {

    print "$gene\t$allele\t$Hash{$gene}{$count}\n";

  }

}

