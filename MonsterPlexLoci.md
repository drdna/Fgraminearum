# _F.graminearum_
Code and data for Population Genomic Studies of _Fusarium graminearum_

## Determine Allele Frequencies of SNPs within the NA1, NA2 and NA3 subgroups of the _F. graminearum_ population
1. Use [FgramAlleleFreqs.pl](/scripts/FgramAlleleFreqs.pl) script to generate a table of allele frequency counts from the FgramHaplotypes.complete.txt file:
```bash
perl FgramAlleleFreqs.pl FgramHaplotypes.complete.txt > FgramAlleleFreqs.txt
```
![FgramAlleleFreqs.png](/data/FgramAlleleFreqs.png)

## BLAST all F. gram genomes against the PH1 reference:
```bash
mkdir FgramBLAST
cd Fgram_unmasked
for f in `ls *nh.fasta`; do blastn -query PH1ref.fasta -subject $f -evalue 1e-20 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' > ../FgramBLAST/PH1.${f/_*/}.BLAST; done
```
## Identify genomic regions that are invariant across all genomes in a given list of strains
Iterate through BLAST alignments of PH1 searched against all other strains and iterate through btop strings to identify blocks with no variants.
```bash
perl PerfectAlignments.pl PH1Ref.fasta FgramNAstrains.txt NewFgramBLAST > FgramInvariant.txt
```
![FgramInvariant.png](/data/FgramInvariant.png)

## Identify chromosome blocks that contain multiple, high frequency SNPs (>=2) within a specified distance of one another (41 nt)
```bash
perl  ClusteredMonsterplexSNPsNewFormat.pl FgramAlleleFreqs.txt 41
```
## Identify high frequency variants that are flanked by invariant priming sites
Grab multiple 100 nt windows that flank each hi-freq SNP site (44 nt central region containing the variants and 28 nt on each end for potential primer annealing sites). Produce Primer3 input file that requests primers with Tm of 54-50, product size range of 800-100 nt and aoids primers that extend into variable sites (positions 29-72).
```bash
perl SelectMonsterplexTargetLoci.pl FgramInvariant.txt FgramAlleleFreqs.txt ClusteredHiFreqSNPs.txt 50 > MPlex_target_candidates.txt
```
## Use primer3 to pick primers (57 < Tm < 60; fragment size range 90-110 bp)
```bash
primer3_core primer3_in.txt > primer3SuggestedPrimers.txt
```
## Grab sequences of primers suggested by primer3 and the sequences of the corresponding loci that would be amplified from PH1 reference genome
```bash
perl PickMPlexPrimers.pl primer3_primer_suggestions.txt Fgram_primer_set1
```
Resulting primer sequence file can be accessed here: [Fgram_primer_set1.csv](/data/Fgram_primer_set1.csv)
Amplified locus sequences can be accessed here: [Fgram_primer_set1.fasta](/data/Fgram_primer_set1.fasta)

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
