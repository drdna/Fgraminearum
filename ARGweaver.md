# Building Ancestral Recombination Graphs (ARGs) for _Fusarium graminearum_
1. Generate ARGweaver sites file using [Generate_ARGsites.pl](/scripts/Generate_ARG.pl) script. We only grabbed data for the first 1 Mb with no data compression.
```bash
Generate_ARGsites.pl FgramARGstrainsPlus.txt FgWardPlusHaplotypes.complete.txt 3 | awk '$1 ~ /NAMES/ || $1 < 1000000' > FgramARGchr3_1-1000000.txt
```
2. Manually edit line 2 in the resulting file to contain info on the region analyzed: REGION 1 Chr3 145118 999584
3. Run ARGweaver on MCC using the [ARGweaver.sh](/scripts/ARGweaver.sh) script:
```bash
sbatch ARGweaver.sh FgramARGchr3_1-1000000.txt 2e-9 1e-11 0
```
4. Extract maximum clade credibility trees from the resulting smc files:
```bash
python3 ARGweaver.py FgramARGchr3_1-1000000.txt_2e-9_1e-11_145118-999584
```
5. Generate tanglegrams using [PlotTanglegrams.py](/scripts/PlotTanglegrams.py) script:
```bash
python3 PlotTanglegrams.py FgramARGchr3_1-1000000.txt_2e-9_1e-11_145118-999584
```

![FG-tanglegram-ML-trees-large.png](/data/FG-tanglegram-ML-trees-large.png)
