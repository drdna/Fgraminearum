# Analysis of Phylogenetic Congrence in Representative Hi-speed versus Lo-speed Genome Regions
## Hi-speed region
1. Retrieve haplotype data for consecutive segments at the start of Chromosome 3:
```bash
for f in $(seq 1 250 5000); do g=$(($f+250)); awk -v var=$f '$3 == "sequence3" {print ">" $2 "_" $1 "\n" substr($4, var, 250)}' FgWardPlusHaplotypes.complete.txt | grep -f NA1NA2NA3strains.txt -A 1 | grep -v ^- > FgramChr3_${f}-${g}.fasta; done
```
2. Build Trees for the each of the selected regions:
```bash
for f in $(seq 1 250 5000); do g=$(($f+250)); /Applications/standard-RAxML-master/raxmlHPC-SSE3 -m BINGAMMA -n FgramChr3_${f}-${g} -s FgramChr3_${f}-${g}.fasta -p 1234 -f a -x 4321 -# 100
```
3. Modify tree filenames for importation into PlotTanglegrams.py script:
```bash
a=0; for f in $(seq 1 250 5000); do g=$(($f+250)); a=$((a+=1)); cp RAxML_bestTree.FgramChr3_${f}-${g} MCCT/MCCT${a}.tre; done
```
4. Add a reference tree in which F. gerlachii, NA1, NA2 and NA3 form discrete clusters:
```bash
mv RAxML_bestree.FgramChr3_27000-32000 MCCT/MCCT0.tre
```
5. Plot Tanglegram:
```bash
python PlotTanglegramsNew.py MCCT/
```
6. This produced the following tanglegram:
![HS-tanglegram-ML-trees.png](/data/HS-tanglegram-ML-trees.png)

## Lo-speed region:
1. Retrieve haplotype data for consecutive segments starting at variant position 21,000 from the start of Chromosome 3. Note: we have to use a larger variant window (2,000) because the NA1,NA2 and NA3 populations have lower SNP density in this region. 
```bash
for f in $(seq 21000 2000 61000); do g=$(($f+2000)); awk -v var=$f '$3 == "sequence3" {print ">" $2 "_" $1 "\n" substr($4, var, 2000)}' FgWardPlusHaplotypes.complete.txt | grep -f NA1NA2NA3strains.txt -A 1 | grep -v ^- > FgramChr3_${f}-${g}.fasta; done
```
2. Build Trees for the each of the selected regions:
```bash
for f in $(seq 21000 2000 61000); do g=$(($f+2000)); /Applications/standard-RAxML-master/raxmlHPC-SSE3 -m BINGAMMA -n FgramChr3_${f}-${g} -s FgramChr3_${f}-${g}.fasta -o Fger_38380,Fger_36905 -p 1234 -f a -x 4321 -# 100
```
3. Modify tree filenames for importation into PlotTanglegrams.py script:
```bash
a=0; for f in $(seq 21000 2000 5000); do g=$(($f+2000)); a=$((a+=1)); cp RAxML_bestTree.FgramChr3_${f}-${g} MCCT/MCCT${a}.tre; done
```
4. Plot Tanglegram:
```bash
python PlotTanglegramsNew.py MCCT/
```
5. This produced the following tanglegram:
![LS-tanglegram-ML-trees.png](/data/LS-tanglegram-ML-trees.png)

### Note how the hi-speed regions fail to resolve the NA1, NA2, NA3, F. gerlachii groupings and, furthermore, show extensive haplotype sharing between all groups.
