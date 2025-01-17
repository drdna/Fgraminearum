# Building Ancestral Recombination Graphs (ARGs) for _Fusarium graminearum_
Assess recombinational history (# events and timing) using ARGweaver.
1. Generate ARGweaver sites file using [Generate_ARGsites.pl](/scripts/Generate_ARG.pl) script. We only grabbed data for the first 1 Mb with no data compression.
```bash
Generate_ARGsites.pl FgramARGstrainsPlus.txt FgWardPlusHaplotypes.complete.txt 3 | awk '$1 ~ /NAMES/ || $1 < 1000000' > FgramARGchr3_1-1000000.txt
```
2. Manually edit line 2 in the resulting file to contain info on the region analyzed: REGION  Chr3 145118 999584 (tab-delimited)
3. Run ARGweaver on MCC using the [ARGweaver.sh](/scripts/ARGweaver.sh) script (ARGweaver.sh \<infile\> \<mutation\> \<recombination\> \<range\>):
```bash
sbatch ARGweaver.sh FgramARGchr3_1-1000000.txt 2e-9 1e-11 1-1000000
```
4. Extract maximum clade credibility trees from the resulting smc files:
```bash
python3 ARGweaver.py FgramARGchr3_1-1000000.txt_2e-9_1e-11_145118-999584
```
5. Add strain names to smc files using SMCnames.pl script:
```bash
for f in `ls *smc.gz`; do perl SMCnames.pl $f; done
```
6. Generate tanglegrams using [PlotTanglegrams_v02.py](/scripts/PlotTanglegrams_v02.py) script:
```bash
python3 PlotTanglegrams.py FgramARGchr3_1-1000000.txt_2e-9_1e-11_145118-999584/
```

![tanglegram-test.png](/data/tanglegram-test.png)

7. Build an Ancestral Recombination Graph
```bash
python SMC2ARG.py FgramARGchr3_1-1000000.txt_2e-9_1e-11_145118-999584/FgramARGchr3_1-1000000.txt_2e-9_1e-11_145118-999584.99.smc
```
8. Plot recombination dates as a function of chromosome position ([RecombinationAge.R](/scripts/RecombinationAge.R)):
![RecombinationAge.png](/data/RecombinationAge.png)

Note that the last 
