# Demonstration of Phylogenetic Ingrongrence in Hi-speed Genome Regions
1. Retrieve haplotype data for consecutive segments as the start of Chromosome 3:
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
![tanglegram-ML-trees.png](/data/tanglegram-ML-trees.png)

