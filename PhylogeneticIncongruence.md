# Analysis of Phylogenetic Congrence in Representative Hi-speed versus Lo-speed Genome Regions
## Hi-speed region
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

#### The first tree is the reference tree which is used to show the consensus grouping for _F. gerlachii_ NA1, NA2 and NA3. Note the extensive and cross-population haplotype sharing between _F. gerlachii_, NA1, NA2 and NA3 as indicated by trees 2 through 20.

8. Next, I built tanglegrams spanning the Hi-speed region at the end of chromosome 3:
![tanglegram-test-1-20.png](/data/tanglegram-test-1-20.png)
![tanglegram-test-21-40.png](/data/tanglegram-test-21-40.png)
![tanglegram-test-41-60.png](/data/tanglegram-test-41-60.png)
![tanglegram-test-61-80.png](/data/tanglegram-test-61-80.png)
![tanglegram-test-81-100.png](/data/tanglegram-test-81-100.png)

## Lo-speed regions:
1. Retrieve haplotype data for consecutive segments starting at variant position 21,000 from the start of Chromosome 3. Note: we have to use a larger variant window (2,000 sites) because the NA1,NA2 and NA3 populations have lower SNP density in this region. 
```bash
for f in $(seq 21000 2000 61000); do g=$(($f+2000)); awk -v var=$f '$3 == "sequence3" {print ">" $2 "_" $1 "\n" substr($4, var, 2000)}' FgWardPlusHaplotypes.complete.txt | grep -f NA1NA2NA3strains.txt -A 1 | grep -v ^- > FgramChr3_${f}-${g}.fasta; done
```
2. Build Trees for the each of the selected regions:
```bash
for f in $(seq 21000 2000 61000); do g=$(($f+2000)); /Applications/standard-RAxML-master/raxmlHPC-SSE3 -m BINGAMMA -n FgramChr3_${f}-${g} -s FgramChr3_${f}-${g}.fasta -o Fger_38380,Fger_36905 -p 1234 -f a -x 4321 -# 100
```
3. Modify tree filenames for importation into PlotTanglegrams.py script:
```bash
a=0; for f in $(seq 21000 2000 61000); do g=$(($f+2000)); a=$((a+=1)); cp RAxML_bestTree.FgramChr3_${f}-${g} MCCT/MCCT${a}.tre; done
```
4. Plot Tanglegram:
```bash
python PlotTanglegramsNew.py MCCT/
```
5. This produced the following tanglegram:
![LD_tanglegram-test.png](/data/LD_tanglegram-test.png)

#### Note how the hi-diversity regions fail to resolve the NA1, NA2, NA3, _F. gerlachii_ groupings and, furthermore, show extensive haplotype sharing between all groups. The NA1 group shows a number of splits across the "lo-diversity" tanglegram because it is an admixed population, whose donor lineages are at varying distances from NA2 and NA3. Otherise, NA1, NA2 and NA3 group members generally cluster together.

Next, I built a tanglegram for one of the lo-speed regions on chromosome 1:
![tanglegram-test-1-24.png](/data/tanglegram-test-1-24.png)

