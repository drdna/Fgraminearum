# _F.graminearum_
Code and data for Population Genomic Studies of _Fusarium graminearum_

## Determine Allele Frequencies of SNPs within the NA1, NA2 and NA3 subgroups of the _F. graminearum_ population
1. Use [FgramAlleleFreqs.pl](/scripts/FgramAlleleFreqs.pl) script to genetrate a table of allele frequency counts from the FgramHaplotypes.complete.txt file:
```bash
perl FgramAlleleFreqs.pl FgramHaplotypes.complete.txt > FgramAlleleFreqs.txt
```
## Identify genomic regions that contain multiple SNPs within a 50 nt target region
Read Allele Frequencies file to identify blocks of SNPs
```bash
perl ClusteredMonsterplexSNPsNewFormat.pl FgramAlleleFreqs.txt 50 > ClusteredHiFreqSNPs.txt
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
## Print out candidate MonsterPlex loci
```bash
perl SelectMonsterplexTargetLoci.pl FgramInvariant.txt FgramAlleleFreqs.txt ClusteredHiFreqSNPs.txt 50 > FgramMPlexTargets.fasta
```
## 
