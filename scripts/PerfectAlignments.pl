#!/usr/bin/perl

die "Usage: perl PerfectAlignments.pl <ref_genome> <strain-list> <blast-dir>\n" if @ARGV != 3;
   
use Data::Dumper;
use FetchGenomeOrig;
 
($refGenome, $strainList, $blastDir) = @ARGV;


#my %StoredStringHash = ();

READ_STRAIN_LIST();

GET_REF_LENGTHS();

WRITE_ALIGN_STRINGS();

READ_DIR();
  
PRINT_ALIGNMENTS();


## SUBROUTINES

sub READ_STRAIN_LIST {
  open(L, $strainList);
  while($L = <L>) {
    chomp($L);
    ($strain, $pop, $include) = split(/\t| +/, $L);
    $StrainHash{$strain} = 1; 
  }
  close L,
}

sub GET_REF_LENGTHS {
  # Get lengths of RefSeq
  
  $GenomeLengthRef = FetchGenomeOrig::getLens($refGenome);
  %GenomeLengthHash = %$GenomeLengthRef;
}    

sub WRITE_ALIGN_STRINGS {
  # write alignment string and reference string
  foreach $Seq (keys %GenomeLengthHash) {
    ($NewSeq = $Seq) =~ s/[^A-Za-z0-9]+$//;	#strip out all non-alphanumerics at end of sequence name
    print "new Seq = $NewSeq\n";	
    $StoredStringHash{$NewSeq} = "1" x $GenomeLengthHash{$Seq};     
  }
}


sub PRINT_STORED_STRINGS {
  # print 1
  @keys = keys %StoredStringHash;
  foreach $key (@keys) {
    print "$key\n$StoredStringHash{$key}\n"
  }

  # print 2
  print "$StoredStringHash{$keys[0]}\n";
}


sub PRINT_ALIGNMENTS {
    foreach $Seq (keys %StoredStringHash) {
        while($StoredStringHash{$Seq} =~ /(1+)/g) {
            $start = $-[0] + 1;
            $end = $+[0];
            print "$Seq\t$start\t$end\n";
        }
    }
}


sub READ_DIR {
  # Read blast directory
  opendir(BLASTDIR, $blastDir) || die "Can't open blast directory for reading\n";
  @Files = readdir(BLASTDIR);
  foreach my $File (@Files) {
    next if $File !~ /BLAST$/;
    ($Q, $S, $suff) = split(/\./, $File);
    if(exists($StrainHash{$S})) {
       ALIGN_STRINGS($File);
    }
  }
}


sub ALIGN_STRINGS {
    my $blastFile = $_[0];
    open(BLAST, "$blastDir/$blastFile") || die "Can't open blast file\n";
    while(my $blast = <BLAST>) {
        chomp($blast);
        @Alignment = split(/\t/, $blast);
        MISALIGNMENTS(\@Alignment);
    }
}     


sub MISALIGNMENTS {
    # mark non-aligned portions of btop string onto alignment string by converting 1s to 0s

    # print "$StoredStringHash{$keys[0]}\n";
    my $AlignRef = @_[0];
    ($qid, $sid, $qs, $qe, $ss, $se, $btop) = @$AlignRef;
    my $pos = $qs;
    my @Align = split(/(\d+)/, $btop);                          # convert BTOP output to list
    shift @Align;
    foreach my $Alignment (@Align) {                            # loop through BTOP list
        if($Alignment =~ /\d+/) {                               # if list member is an integer, increment position 
            $pos += $Alignment
        }
        elsif($Alignment =~ /\w+/) {                            # if list members are alphabetic, analyze the mismatch and increment reference position if non-gap
            @mismatches = split(//, $Alignment);
            for($i=0; $i <= @mismatches-2; $i+=2) {
                if($mismatches[$i] ne '-') {
                    substr($StoredStringHash{$qid}, $pos-1, 1) =~ tr/1/0/;
                    $pos ++
                }
            }
        }
    }
}

