# Masselink2024_TailRegen

This repository contains the scripts necessary to reproduce the clone calling published in **Somite-independent regeneration of the axolotl primary body axis** (https://doi.org/10.1101/2024.01.31.577464)



The scripts are numbered in the order the should be exectuted:
  1. 01_ViralBarcodes_scRNAseqQuantification.sh (should be executed on shell
  2. 02_scRNAseqProcessing_Rep1.r (should be exectued on R)
  3. 03_scRNAseqProcessing_Rep2.r (should be exectued on R)
  4. 04_BarcodeExtraction.nf (executed on nextflow per replicate fq.gz file)
  5. 05_CloneCalling.r (executed on R)

## Additional files 

1. 99_SessionInfo.txt - File containing the sessionInfo() of the R session this analysis was performed on
 2. Dockerfile - Docker recipe to generate the container necesary for the nextflow script (needs to be added to the script separately)

## Important notes

- The version of SeuratDisk used was 0.0.0.9019 (not included in the sessionInfo file as this was done separately)
- The version of kb this datasets were quantified was 0.24.4
  