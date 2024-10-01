use strict;
use warnings;

# Open the input file
my $filename = $ARGV[0];
open my $fh, '<', $filename or die "Could not open file '$filename': $!";

# Prepare to store column indices and zero counts
my ($numlines, %Positions, %Counts);
    
while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^sequence/;

    if ($line =~ /^sites/) {
        my @positions = split ' ', $line;
        shift @positions;
        my $seqid = shift @positions;
        @{$Positions{$seqid}} = @positions;
    }
    else {
        my ($sampleid, $popinfo, $seqid, $haplotype) = split ' ', $line;
        $numlines++;  # Increment line count

        for (my $i = 0; $i < length($haplotype); $i++) {
            if (substr($haplotype, $i, 1) eq '0') {
                my $position = $Positions{$seqid}[$i];
                $Counts{$seqid}{$position}++;
            }
        }
    }
}

# Determine frequency in each column
my $numseqs = keys %Counts;
foreach my $seqid (sort {$a cmp $b} keys %Counts) {
    foreach my $position (sort {$a <=> $b} keys %{$Counts{$seqid}}) {
        my $freq = $Counts{$seqid}{$position} / ($numlines / $numseqs);
        printf("%s\t%d\t%.2f\n", $seqid, $position, $freq);
    }
}

# Close the file
close $fh;

