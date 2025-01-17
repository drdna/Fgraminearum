#!/usr/bin/perl

## Add subroutine to count N's in columns and filter columns if there are too many.

## Easily done by splicing @Genotype array

## Just iterate through the outfiles twice: count first, filter & print second

##############################
#
# Generate_ARGsites.pl
#
# written by Mark L. Farman
#
# Purpose: Read SNPcaller outfiles, score all SNP positions, and record haplotypes for all individuals in a specified list
#
# Note: Must be re-run every time a new strain is added (unavoidable)
#
##############################

#use strict;

#use warnings;

use MasterList;

my ($strains, $haplFile, $Chr) = @ARGV;

$StrainsHashRef = MasterList::strains($strains);

%StrainsHash = %{$StrainsHashRef};

open(H, $haplFile) || die "Can't open haplotype file\n";

while($H = <H>) {

  chomp($H);
  
  @H = split(/\t/, $H);

  if($H[0] =~ /sites/ && $H[1] =~ /$Chr$/) {
    @SitesList = split(/ /, $H[2]);
  }
  elsif($H[0] =~ /snps/ && $H[1] =~ /$Chr$/) {
    @{$SNPsHash{$H[1]}} = split(/ /, $H[2]);
    $Ref = join '', map {substr($_, 0, 1)} @{$SNPsHash{$H[1]}};
    $Alt = join '', map {substr($_, 1, 1)} @{$SNPsHash{$H[1]}};
  }
  elsif ($H[2] =~ /$Chr$/ && exists($StrainsHash{$H[0]})) {
    push @Names, $H[0];
    while($H[3] =~ /(1)/g) {
      substr($H[3], $-[0], 1, substr($Ref, $-[0], 1))
    }
    while($H[3] =~ /(0)/g) {
      substr($H[3], $-[0], 1, substr($Alt, $-[0], 1))
    }
    @Hapl = split(//, $H[3]);
    push @Haplotypes, [@Hapl];
  }
}

if($ARGV[3] =~ /fasta/) {
  for $i (0..@Names-1) {
    print ">$Names[$i]\n";
    print join ('', @{$Haplotypes[$i]}), "\n";
  }
  exit
}  

print join ("\t","NAMES", @Names), "\n";

for $i (0..@{$Haplotypes[0]}-1) {
  %NuclHash = ();
  my @Line= map {$_->[$i]} @Haplotypes;
  my $Line = join ('', @Line);
  while($Line =~ /([AGTC])/g) {
    $NuclHash{$1} = 1;
  }
  $numNucl = scalar keys %NuclHash;
  if($numNucl > 1) {
    print "$SitesList[$i]\t$Line\n";
  }
} 
 
    
  

