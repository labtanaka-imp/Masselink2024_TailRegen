
#Generate Index for quantification
kb ref -i AmexG_v6.DD \
-g AmexG_v6.DD.T2G \
--workflow lamanno \
--overwrite \
-f1 AmexG_v6.DD.cDNA.fa \
-f2 AmexG_v6.DD.intron.fa \
-c1 AmexG_v6.DD.cDNA.T2C \
-c2 AmexG_v6.DD.intron.T2C \
AmexG_v6.DD.corrected.round2.FULL.mtMasked.mtIncluded.fa \
AmexT_v47.FULL.IDMod.gtf

#Quantify Replicate1
kb count -i ./AmexG_v6.DD \
-g ./AmexG_v6.DD.T2G \
-x 10XV2 \
-o 197397 \
--workflow lamanno \
--overwrite \
--mm \
--filter bustools \
-t 32 \
-m 390G \
--h5ad \
-c1 ./AmexG_v6.DD.cDNA.T2C \
-c2 ./AmexG_v6.DD.intron.T2C \
./197397_S1_L004_R1_001.fastq.gz \
./197397_S1_L004_R2_001.fastq.gz

#Quantify Replicate2
kb count -i ./AmexG_v6.DD \
-g ./AmexG_v6.DD.T2G \
-x 10XV2 \
-o 208683 \
--workflow lamanno \
--overwrite \
--mm \
--filter bustools \
-t 32 \
-m 390G \
--h5ad \
-c1 ./AmexG_v6.DD.cDNA.T2C \
-c2 ./AmexG_v6.DD.intron.T2C \
./208683_S1_L003_R1_001.fastq.gz \
./208683_S1_L003_R2_001.fastq.gz

#Convert the files to seurat object
