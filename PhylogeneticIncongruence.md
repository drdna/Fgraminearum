# Analysis of Phylogenetic Congrence in Representative Hi-speed versus Lo-speed Genome Regions
## Hi-diversity region
1. Retrieve haplotype data for consecutive segments at the start of Chromosome 3:
```bash
for f in $(seq 1 250 5000); do g=$(($f+250)); awk -v var=$f '$3 == "sequence3" {print ">" $2 "_" $1 "\n" substr($4, var, 250)}' FgWardPlusHaplotypes.complete.txt | grep -f NA1NA2NA3strains.txt -A 1 | grep -v ^- > FgramChr3_${f}-${g}.fasta; done
```
2. Build Trees for the each of the selected regions, rooted on _F. gerlachii_:
```bash
for f in $(seq 1 250 5000); do g=$(($f+250)); /Applications/standard-RAxML-master/raxmlHPC-SSE3 -m BINGAMMA -n FgramChr3_${f}-${g} -s FgramChr3_${f}-${g}.fasta -p 1234 -f a -x 4321 -# 100
```
3. Modify tree filenames for importation into PlotTanglegrams.py script:
```bash
a=0; for f in $(seq 1 250 5000); do g=$(($f+250)); a=$((a+=1)); cp RAxML_bestTree.FgramChr3_${f}-${g} MCCT/MCCT${a}.tre; done
```
4. Add a reference tree in which _F. gerlachii_, NA1, NA2 and NA3 form discrete clusters:
```bash
mv RAxML_bestree.FgramChr3_27000-32000 MCCT/MCCT0.tre
```
5. Plot Tanglegram:
```bash
python PlotTanglegramsNew.py MCCT/
```
6. This produced the following tanglegram:
![HD_tanglegram-test.png](/data/HD_tanglegram-test.png)

#### The first tree is the reference tree which is used to show the consensus grouping for _F. gerlachii_ NA1, NA2 and NA3. Note the extensive, cross-population haplotype sharing between _F. gerlachii_, NA1, NA2 and NA3 as indicated by trees 2 through 20.

8. Next, I built tanglegrams spanning the Hi-speed region at the end of chromosome 3:
![tanglegram-test-1-20.png](/data/tanglegram-test-1-20.png)
![tanglegram-test-21-40.png](/data/tanglegram-test-21-40.png)
![tanglegram-test-41-60.png](/data/tanglegram-test-41-60.png)
![tanglegram-test-61-80.png](/data/tanglegram-test-61-80.png)
![tanglegram-test-81-100.png](/data/tanglegram-test-81-100.png)

#### Again, note the poor phylogenetic congruence and the failure of _F. gerlachii_, NA1, NA2 and NA3 to form reliable groupings. Note also there are very few regions where the consensus phylogenetic relationships exhibited by the lo-speed regions are recovered.

## Lo-diversity regions:
1. Retrieve haplotype data for consecutive segments starting at variant position 21,000 from the start of Chromosome 3. Note: we have to use a larger variant window (2,500 sites) because the NA1,NA2 and NA3 populations have lower SNP density in this region. 
```bash
for f in $(seq 21000 2500 61000); do g=$(($f+2500)); awk -v var=$f '$3 == "sequence3" {print ">" $2 "_" $1 "\n" substr($4, var, 2500)}' FgWardPlusHaplotypes.complete.txt | grep -f NA1NA2NA3strains.txt -A 1 | grep -v ^- > FgramChr3_${f}-${g}.fasta; done
```
2. Build Trees for the each of the selected regions:
```bash
for f in $(seq 21000 2500 61000); do g=$(($f+2500)); /Applications/standard-RAxML-master/raxmlHPC-SSE3 -m BINGAMMA -n FgramChr3_${f}-${g} -s FgramChr3_${f}-${g}.fasta -o Fger_38380,Fger_36905 -p 1234 -f a -x 4321 -# 100
```
3. Modify tree filenames for importation into PlotTanglegrams.py script:
```bash
a=0; for f in $(seq 21000 2500 61000); do g=$(($f+2500)); a=$((a+=1)); cp RAxML_bestTree.FgramChr3_${f}-${g} MCCT/MCCT${a}.tre; done
```
4. Plot Tanglegram:
```bash
python PlotTanglegramsNew.py MCCT/
```
5. This produced the following tanglegrams:
![tanglegram-test-1-16.png](/data/tanglegram-test-1-16.png)
![tanglegram-test-17-32.png](/data/tanglegram-test-17-32.png)

#### Note how the lo-diversity regions exhibit good phylogenetic congruence. In particular, NA3 is consistently closest to _F gerlachii_ and shows no phylogenetic conflict with NA1 or NA2. The NA1 group shows a number of splits across the "lo-diversity" tanglegram because it is an admixed population, whose donor lineages are at varying distances from NA2 and NA3. Otherise, NA1, NA2 and NA3 group members generally cluster together, with the rare exception of recombinant individuals (n=5).

Next, I built a tanglegram for one of the lo-speed regions on chromosome 1:
![tanglegram-test-1-24.png](/data/tanglegram-test-1-24.png)

#### This also exhibited good phylogenetic congruence and, again, pointed to the existence of two discrete haplotypes at several loci within the NA1 population.
