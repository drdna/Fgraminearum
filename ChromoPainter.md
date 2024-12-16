# ChromoPainting of foreign DNA introgressions
1. Use GenerateCP.pl script to generate input files from FgWardPlusHaplotypes.complete.txt file:
```bash
perl GenerateCP.pl FgWardPlusHaplotypes.complete.txt sequence3 FgWardPlusChr3.cp
```
2. Edit [FgWardPlusChr3.pop](/data/FgWardPlusChr3.pop) file to specify donor and recipient populations.
3. Run ChromoPainterV2:
```bash
ChromoPainterV2 -j -g FgWardPlusChr3.cp -r FgWardPlusChr3.recomb -t FgWardPlusChr3.idfile -f FgWardPlusChr3.pop 0 0 -o FgWardPlusChr3all_out -i 10 -ip -b
```
4. Subsample the dataset using the [ChromoPaint_to_R_compressed.pl](/scripts/ChromoPaint_to_R_compressed.pl) script:
```bash
perl ChromoPaint_to_R_compressed.pl FgWardPlusChr3all_out.copyprobsperlocus.out.gz
```
5. Copy the data into a dedicated directory:
```bash
mkdir FgramCopyProbs
mv *copyprobsperlocus.out.gz FgramCopyProbs
```
6. Plot foreign introgressions using the [ForeignIntrogressions.R](/scripts/ForeignIntrogressions.R) script:
