#!/usr/bin/perl

die "Usage: perl ClusteredMonsterplexSNPsNewFormat.pl <allele-frequencies-list> <cluster-size>\n" if @ARGV < 2;

use strict;
use warnings;

my ($infile, $clustersize) = @ARGV;

open my $fh, '<', $infile or die "Could not open file: $!";

my %sequences;

# Read the file and store positions by sequence
while (my $line = <$fh>) {
    chomp $line;
    my ($seq, $pos, $freq) = split '\t| +', $line;  # Split line into sequence, position and frequency
#   print "$seq, $pos, $freq\n";
    next unless $freq > 0.1 && $freq < 0.9;
    push @{ $sequences{$seq} }, $pos;    # Push positions into an array for each sequence
}

close $fh;

# Iterate over sequences and find clusters
foreach my $seq (keys %sequences) {
    my @positions = sort { $a <=> $b } @{ $sequences{$seq} };  # Sort positions
    my @clustered;

    for (my $i = 0; $i < @positions; $i++) {
        if (!@clustered) {
            push @clustered, $positions[$i];  # Start a new cluster
        } else {
            if ($positions[$i] - $clustered[-1] <= $clustersize) {
                push @clustered, $positions[$i];  # Continue the cluster
            } else {
                # Print the start and end positions of current cluster and start a new one
#                print "$seq\t$clustered[0]\t$clustered[-1]\n" if @clustered > 1;
                print join(" ", $seq, @clustered), "\n\n" if @clustered > 1;	# print whole array of sites
                @clustered = ($positions[$i]);  # Start a new cluster
            }
        }
    }
    # Print any remaining clustered positions
    if (@clustered) {
        print "$seq\t@clustered\n" if @clustered > 1;
#        print "$seq\t$clustered[0]\t$clustered[-1]\n" if @clustered > 1;
    }
}

