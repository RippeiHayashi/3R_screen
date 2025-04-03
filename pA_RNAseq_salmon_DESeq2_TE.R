### the code was used to perform differential gene expression analysis using salmon quant.sf files as inputs

pdf.options(useDingbats = FALSE)

library("readr")
library("tximport")
library("tximportData")
library("DESeq2")
library("ggplot2")
library("ggrepel")
library("EnhancedVolcano")

setwd("<working_directory>")
TMP="<working_directory>"
COMPARISON="SpnE_pA_RNAseq"
CONTRAST="genotype"
### when contrasting ZF and WT
CONTRAST_A="ZF"
CONTRAST_B="WT"

### input files
v=c(TMP,"/",COMPARISON,"/",COMPARISON,"_salmon_quant_files_TE.txt")
SALMON.FILES=paste(v,collapse="")
v=c(TMP,"/",COMPARISON,"/",COMPARISON,"_library_names.txt")
LIBS=paste(v,collapse="")
v=c(TMP,"/",COMPARISON,"/",COMPARISON,"_library_categories.txt")
CATEGORIES=paste(v,collapse="")

### output files
v=c(TMP,"/",COMPARISON,"/",COMPARISON,"_normalised_counts_TE.csv")
normalised_counts=paste(v,collapse="")
v=c(TMP,"/",COMPARISON,"/",COMPARISON,"_",CONTRAST_A,"_over_",CONTRAST_B,"_differential_gene_expression_TE.txt")
DE_TE=paste(v,collapse="")
v=c(TMP,"/",COMPARISON,"/",COMPARISON,"_",CONTRAST_A,"_over_",CONTRAST_B,"_plotMA_TE.labeled.pdf")
plotMA_name_TE=paste(v,collapse="")
mains=c("MA plot, ",COMPARISON,", NonSig:gray&black, Sig(q-value<0.05):pink&red, TE:black&red")
mains_name=paste(mains,collapse="")

### DESeq2 set up and output results table --- START ---
### import a file with the path to quant.sf file per line
salmon_files <- read.table(file=SALMON.FILES,header=FALSE)
### convert a data.frame column to a vector
salmon_files_vector <- as.vector(salmon_files[,1])
### import a file with library name per line
file_names <- read.table(file=LIBS,header=FALSE)
### put filenames
names(salmon_files_vector) <- file_names[,1]

### category file
sample_categories <- read.table(file=CATEGORIES, header=TRUE)
rownames(sample_categories) <- file_names[,1]

### file that matches transcript ID to gene ID
tx2gene <- read_csv("fbgn_fbtr_2019_03_plus_122TE.modified.txt")

### measure gene_based read counts
tx <- tximport(salmon_files_vector, type="salmon",tx2gene=tx2gene)

### carry out DESeq
ddsTxi <- DESeqDataSetFromTximport(tx,colData = sample_categories,design = ~ genotype)
ddsTxi <- DESeq(ddsTxi)

### normalised counts
foo <- counts(ddsTxi, normalized = TRUE)
colSums(foo)
write.csv(foo,file=normalised_counts)
### discard lowly expressed genes
CUTOFF=10
CUTOFF=as.numeric(CUTOFF)
keep <- rowSums(foo) >= CUTOFF
### carry out DESeq
ddsTxi.f <- ddsTxi[keep,]

### perform DESeq2::results
res.f <- results(ddsTxi.f, contrast=c(CONTRAST, CONTRAST_A, CONTRAST_B))
write.table(res.f,file=DE_TE)

### take TE data
foo.f <- counts(ddsTxi.f, normalized = TRUE)
keep_TE <- grep("TE", rownames(foo.f), value = FALSE)
ddsTxi.TE <- ddsTxi.f[keep_TE,]
res.TE <- results(ddsTxi.TE, contrast=c(CONTRAST, CONTRAST_A, CONTRAST_B))

### select TEs for labeling
topTE <- c("TE:HeT-A", "TE:Burdock", "TE:TAHRE", "TE:gypsy12", "TE:GATE", "TE:Max-element", "TE:invader2")

### make MA plot
mainy=c("log2 fold change in ",CONTRAST_A," over ",CONTRAST_B)
mainy_name=paste(mainy,collapse="")
pdf(file=plotMA_name_TE,width=10,height=6)
plotMA(res.f, ylim=c(-8,8), xlim=c(20,2000000),cex = 0.5, alpha=0.05, colNonSig = "lightgray", colSig = "lightpink", ylab="")
par(new=T)
plotMA(res.TE, ylim=c(-8,8), xlim=c(20,2000000),cex = 1, alpha=0.05, colNonSig = "black", colSig = "red", ylab=mainy_name, main=mains_name)
with(res.TE[topTE, ], {
         points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
         text(baseMean, log2FoldChange, topTE, pos = 2, col = "dodgerblue")
})
dev.off()

sessionInfo()

R version 4.3.1 (2023-06-16)
attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] EnhancedVolcano_1.20.0      ggrepel_0.9.4               ggplot2_3.4.4               DESeq2_1.42.0              
 [5] SummarizedExperiment_1.32.0 Biobase_2.62.0              MatrixGenerics_1.14.0       matrixStats_1.1.0          
 [9] GenomicRanges_1.54.1        GenomeInfoDb_1.38.1         IRanges_2.36.0              S4Vectors_0.40.1           
[13] BiocGenerics_0.48.1         tximportData_1.30.0         tximport_1.30.0             readr_2.1.4                

loaded via a namespace (and not attached):
 [1] generics_0.1.3          utf8_1.2.4              SparseArray_1.2.2       bitops_1.0-7            lattice_0.22-5         
 [6] hms_1.1.3               magrittr_2.0.3          grid_4.3.1              Matrix_1.6-1.1          fansi_1.0.5            
[11] scales_1.2.1            codetools_0.2-19        abind_1.4-5             cli_3.6.1               rlang_1.1.1            
[16] crayon_1.5.2            XVector_0.42.0          munsell_0.5.0           withr_2.5.2             DelayedArray_0.28.0    
[21] S4Arrays_1.2.0          tools_4.3.1             parallel_4.3.1          tzdb_0.4.0              BiocParallel_1.36.0    
[26] dplyr_1.1.3             colorspace_2.1-0        locfit_1.5-9.8          GenomeInfoDbData_1.2.11 vctrs_0.6.4            
[31] R6_2.5.1                lifecycle_1.0.3         zlibbioc_1.48.0         pkgconfig_2.0.3         pillar_1.9.0           
[36] gtable_0.3.4            glue_1.6.2              Rcpp_1.0.11             tidyselect_1.2.0        tibble_3.2.1           
[41] rstudioapi_0.16.0       compiler_4.3.1          RCurl_1.98-1.13        
