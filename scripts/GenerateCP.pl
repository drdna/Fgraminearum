#!/usr/bin/perl


my ($filename, $chromosome, $cpfile) = @ARGV;

die "CP outfile must end in .cp\n" if $cpfile !~ /cp$/;

my $lineCount = `wc -l < $filename`;
chomp($lineCount);
$lineCount =~ s/^\s+//;
my $seqCount = `grep -c sites < $filename`;
chomp($seqCount);
my $numStrains = ($lineCount - (2 * $seqCount))/$seqCount;

open(CP, '>', $cpfile) || die "Can't open outfile\n";
print CP "$numStrains\n";

($strainIdFile = $cpfile) =~ s/\.cp/.idfile/;
open(ID, '>', $strainIdFile) || die "Can't open strain id file\n";

($recombfile = $cpfile) =~ s/\.cp/.recomb/;
open(RECOMB, '>', $recombfile) || die "Can't open recombfile\n";
print RECOMB "start.pos recom.rate.perbp\n";


open(H, $filename);
while($H = <H>) {
  chomp($H);
  @H = split(/\t| +/, $H);
  $len = @H;
  next if $H[0] =~ /^snps/;
  if($H[0] =~ /^sites/) {
    my $sequence = $H[1];
    if($sequence =~ /$chromosome$/) {
      $numSites = $len - 2;
      print CP "$numSites\n";
      print CP join (" ", "P", @H[2..$#H]), "\n";
      foreach $value (@H[2..$#H]) {
        $count ++;
        if($count == 1) {
          print RECOMB "$value" || die "$!\n"
        }
        else {
          print RECOMB " 0.0000007\n$value"
        } 
      }  
      print RECOMB " 0\n";
    }
  }
  elsif($H[2] =~ /$chromosome$/) {
    $popHash{$H[1]} = 1;
    print CP "$H[3]\n";
    print ID "$H[0] $H[1] 1\n"
  }
}

close CP;
close RECOMB;
close H;
close ID;

($popfile = $cpfile) =~ s/\.cp/.pop/;
open(POP, '>', $popfile) || die "Can't open popfile\n";
foreach my $pop (keys %popHash) {
  print POP "$pop\n";
}
close POP;
