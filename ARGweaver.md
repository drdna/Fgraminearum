# Building Ancestral Recombination Graphs (ARGs) for _Fusarium graminearum_
Assess recombinational history (# events and timing) using ARGweaver.
1. Generate ARGweaver sites file using [Generate_ARGsites.pl](/scripts/Generate_ARGsites.pl) script. We only grabbed data for the first 1 Mb with no data compression.
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

Note: the last six trees are in the lo-diversity region of the chromosome.

7. Plot leaf traces
```bash
arg-layout FgramARGchr3_1-1000000.txt_2e-9_1e-11_145118-999584/FgramARGchr3_1-1000000.txt_2e-9_1e-11_145118-999584.199.layout FgramARGchr3_1-1000000.txt_2e-9_1e-11_145118-999584/FgramARGchr3_1-1000000.txt_2e-9_1e-11_145118-999584.199.smc
```
8. Load layout into [PlotLeafTraces.R](/script/PlotLeafTraces.R) script:

![FgramARGchr3LeafTraces.png](/data/FgramARGchr3LeafTraces.png)

9. Build an Ancestral Recombination Graph
```bash
python SMC2ARG.py FgramARGchr3_1-1000000.txt_2e-9_1e-11_145118-999584/FgramARGchr3_1-1000000.txt_2e-9_1e-11_145118-999584.99.smc
```
10. Plot recombination dates as a function of chromosome position ([RecombinationAge.R](/scripts/RecombinationAge.R)):

![RecombinationAge.png](/data/RecombinationAge.png)

Note: the last six trees are in the lo-diversity region of the chromosome.

11. Extend the ARGs further into the lo-diversity region (up to position 3,000,000):

The portion from 1,000,000 to 2,000,000 produced the following tanglegram, which reveals a hi-diversity/hi-recombination region spanning positions 1,700,000 to 1,850,000: 

![LD_tanglegram-test2.png](/data/LD_tanglegram-test2.png)

and the entire block (1,000,000 to 3,000,000) produced the following distribution of recombination dates:
 
![RecombinationAge2.png](/data/RecombinationAge2.png)

