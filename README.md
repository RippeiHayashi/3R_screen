# 3R_screen
Code for analysing the poly-A selected RNA sequencing and small RNA sequencing data from Drosophila ovaries.

# 3R_screen_RNAseq_analysis.sh
Code was used to analyse the poly-A selected paired end RNA sequencing data from Drosophila ovaries. salmon was used to map reads to the Drosophila transcriptome including transposon consensus sequences (output file: quant.sf). 

# 3R_screen_sRNAseq_analysis.sh
Code was used to analyse the small RNA sequencing data from Drosophila ovaries.

# pA_RNAseq_salmon_DESeq2_TE.R
Code was used to run the differential gene expression analysis using DESeq2, taking the transcript count files from salmon.

# pingpong_linkage.R
Code to measure the linkage score of the ping-pong piRNA biogenesis.

# plotting_TE-piRNAs.R
Code to plot the 5' ends of piRNAs mapping to transposon sequences.

# 122_dmel_TE_Senti2015.fa
Consensus sequences of the 122 Drosophila transposons, taken from Senti, G&D, 2015 (PMID: 26302790).

# fbgn_fbtr_2019_03_plus_122TE.modified.txt
A table for converting the transcript ID from FlyBase (dm6 r6.31) to the gene ID, including transpospon entries.
