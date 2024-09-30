# Fgraminearum
Code and data for Population Genomic Studies of Fusarium graminearum

## Determining Allele Frequencies of SNPs within _F. graminearum_ population
1. Read haplotypes file and format SNPs into columns:
```bash
awk 'BEGIN {OFS=" "} 
{
    if ($1 ~ /sites/) {
        for (i=2; i<=22; i++) {
            printf "%s%s", $i, (i < 22 ? OFS : "\n")
        }
    } else {
        printf "%s ", $3;  # Print the third field followed by a space
        for (j=1; j<=20; j++) {
            printf "%s%s", substr($4, j, 1), (j < 20 ? OFS : "\n")
        }
    }
}' FgramHaplotypes.complete.txt > sequences.txt
```
