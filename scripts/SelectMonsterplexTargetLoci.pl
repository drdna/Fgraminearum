#!/usr/bin/perl
#use strict;
#use warnings;
use FetchGenomeOrig;

die "Usage: SelectMonsterplexTargetLoci.pl <non-polymorphic-alignStrings> <clustered-snp-sites> <ref-genome>\n" if @ARGV != 3;


($NPalignString, $clusteredFrequentSNPs, $RefGenome) = @ARGV;

# Open the input file for reading (replace 'input.txt' with your file name)
open my $fh, '<', $NPalignString or die "Cannot open file\n";


GET_GENOME();

FIND_INVARIANT_BLOCKS();

FIND_TARGETs();


sub GET_GENOME {

  $genomeHashRef = FetchGenomeOrig::getSeqs($RefGenome);
  %genomeHash = %$genomeHashRef;
  @keys = keys %genomeHash;
#  print "$genomeHash{unplaced}\n";

}

sub FIND_INVARIANT_BLOCKS {

# Read the Invariant sequence blocks file, line by line
# Write invariant sites as 1s on a string of 0s
# e.g. 111110111011100111010111110111111111110111111111011011111111111111

    while (my $line = <$fh>) {
        chomp $line;  # Remove newline characters
        my ($chr, $start, $end) = split /\t| +/, $line;  # Split the line into chromosome, start and end of invariant block
        next if $chr !~ /Chr/;
        $Counts{$chr} = "0" unless exists($Counts{$chr});	
        $Counts{$chr} .= "0" x ($end - length($Counts{$chr})) if ($end > length($Counts{$chr}));
        substr($Counts{$chr}, $start -1, 1 + $end - $start) =~ tr/0/1/;            
    }
    # Close the file handle
    close $fh;
}


sub FIND_TARGETs {
# iterate through clustered SNPs file and find loci that are flanked by invariant (primer-binding) regions
# target PCR product is 100 bp (+/- 10)
# need 35 nt invariant at each end

    open(CS, "$clusteredFrequentSNPs") || die "Can't open file: $!\n";
    while($CS = <CS>) {
        chomp($CS);
        $lastChr = $Chr;
        ($Chr, $site1, $site2) = split(/\t/, $CS);
        $windowEnd = 0 if $Chr ne $lastChr;
#        next if $site1 < $windowEnd;
        $Chr =~ s/sequence/Chr/;
        $midpoint = $site1 + int($site2 - $site1 + 0.5);
        ($windowStart, $windowEnd) = ($midpoint - 50, $midpoint + 49);
        $primerRegion1 = substr($Counts{$Chr}, $windowStart, 28);
        if($primerRegion1 !~ /0/) {
             $primerRegion2 = substr($Counts{$Chr}, $windowEnd-28-1, 28);
             if($primerRegion2 !~ /0/) {
                 $primerLocus = substr($genomeHash{$Chr}, $windowStart, 100);
                 print "SEQUENCE_ID=$Chr"."_"."$site1"."_"."$site2\nSEQUENCE_TEMPLATE=$primerLocus\n". 
                 "PRIMER_MIN_TM=54\nPRIMER_MAX_TM=60\nPRIMER_PRODUCT_SIZE=80-100\nPRIMER_NUM_RETURN=1\n".
                 "SEQUENCE_INTERNAL_EXCLUDED_REGION=29,44\n=\n"; 
             }
        }
    }
}

