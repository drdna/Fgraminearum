# _F.graminearum_
Code and data for Population Genomic Studies of _Fusarium graminearum_

## Determine Allele Frequencies of SNPs within the NA1, NA2 and NA3 subgroups of the _F. graminearum_ population
1. Use [FgramAlleleFreqs.pl](/scripts/FgramAlleleFreqs.pl) script to genetrate a table of allele frequency counts from the FgramHaplotypes.complete.txt file:
```bash
perl FgramAlleleFreqs.pl FgramHaplotypes.complete.txt > FgramAlleleFreqs.txt
```
## BLAST all F. gram genomes against the PH1 reference:
```bash
mkdir FgramBLAST
cd Fgram_unmasked
for f in `ls *nh.fasta`; do blastn -query PH1ref.fasta -subject $f -evalue 1e-20 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' > ../FgramBLAST/PH1.${f/_*/}.BLAST; done
```
## Identify genomic regions that are invariant across all genomes in a given list of strains
```bash
perl PerfectAlignments.pl PH1Ref.fasta FgramNAstrains.txt NewFgramBLAST > FgramInvariant.txt
```
## Identify high frequency variants that are flanked by invariant priming sites
Grab multiple 91 nt windows that flank each hi-freq SNP site (41 nt central region containing the variant and 25 nt on each end for potential primer annealing sits). Slide windows across the locus so that variant site is always >= 35 nt from window/start end.
```bash
perl SelectMonsterplexTargetLoci.pl FgramInvariant.txt FgramAlleleFreqs.txt ClusteredHiFreqSNPs.txt 50 > MPlex_target_candidates.txt
```
## Use primer3 to pick primers (57 < Tm < 60; fragment size range 81-91 bp)
```bash
awk '{print "SEQUENCE_ID=" $1 "\nSEQUENCE_TEMPLATE=" $2 "\nPRIMER_MIN_TM=57\nPRIMER_MAX_TM=60\nPRIMER_PRODUCT_SIZE=81-91\nPRIMER_NUM_RETURN=1\n="}' MPlex_target_candidates.txt > primer3_in.txt
primer3_core primer3_in.txt > primer3_primer_suggestions.txt
```
## Select primers based on the sequences most frequently picked by primer3
```bash
perl PickMPlexPrimers.pl primer3_primer_suggestions.txt > Picked_primers.txt
```
Resulting file can be accessed here: [Picked_primers.txt](/data/Picked_primers.txt)

## Plot primer distribution
Use [PlotPrimerSites.R](/scripts/PlotPrimerSites.R) script to plot primer sites on chromosomes.

![MonsterPlexTargets.png](/data/MonsterPlexTargets.png)

## Develop MonsterPlex Primers for current multi-locus genotyping markers:
1. Blast existing markers against the NA1, NA2 and N3 population members:
```bash
for f in `ls MLG_FASTA/*fasta`; do while read col1 col2 col3; do blastn -query $f -subject Documents/FGRAM/FgramFasta/${col1}_nh_masked.fasta -outfmt '6 qseqid sseqid qlen qstart qend sstart send btop' | awk '$5 - $4 == $3 - 1 {print $1, $8}'; done < FgramNAstrains.txt; done > MLGallelesAll.txt
```
2. Use [SortMLGAlleles.pl](/scripts/SortMLGAlleles.pl) script to determine frequency of each allelic variant:
```bash
perl SortMLGAlleles.pl MLGallelesAll.txt  | sort -k1,1 -k3n > MLGalleleCounts.txt
```
Resulting file can be accessed at [MLGalleleCounts.txt](/data/MLGalleleCounts.txt)
3. Design primers to amplify the most frequent polymorphic sites (target = two per locus) 
