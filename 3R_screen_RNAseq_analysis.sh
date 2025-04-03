### the code was used to analyse the poly-A selected RNA sequencing data from Drosophila ovaries.

### software used
### fastx_toolkit
http://hannonlab.cshl.edu/fastx_toolkit/
### fastq_pair
https://github.com/linsalrob/fastq-pair
### salmon/1.1.0
https://github.com/COMBINE-lab/salmon/releases

### define the variables
references="<directory where the genome sequence file and the gtf file are kept>"
raw_data="<directory where the fastq files are kept>"
analysis="<directory where the output files are kept>"
lib_pA="<name of the poly-A+ RNA seq library>"
mkdir -p ${analysis}/${lib_pA}

### below are the poly-A selected RNA sequencing libraries processed in this script:
w1118_ovaries_pA_rep1
w1118_ovaries_pA_rep2
w1118_ovaries_pA_rep3
spnE_ZF_mutant_ovaries_pA_rep1
spnE_ZF_mutant_ovaries_pA_rep2
spnE_ZF_mutant_ovaries_pA_rep3
spnE_null_mutant_ovaries_pA_rep1
spnE_null_mutant_ovaries_pA_rep2
spnE_null_mutant_ovaries_pA_rep3

### structure of RNA seq libraries
5’ AATGATACGGCGACCACCGAGATCTACAC[TCTTTCCCTACACGACGCTCTTCCGATCT]--- [R1 ->] --- [<- R2] ---[AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC] ---NNNNNN--- ATCTCGTATGCCGTCTTCTGCTTG 3’
R1: reverse complement
R2: sense


### step 0: making a salmon index containing dm6 r6.31 transcriptome and 122 TE sequences from Senti et al., 2015 (also made available in the git repository)
mkdir -p ${references}/dm6/indices
wget --directory-prefix="${references}/dm6" https://ftp.flybase.net/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-transcript-r6.31.fasta.gz
zcat ${references}/dm6/dmel-all-transcript-r6.31.fasta.gz > ${references}/dm6/dmel-all-transcript-r6.31.fasta
### build index for transcriptome + 122TE
salmon index -t <(cat ${references}/dm6/dmel-all-transcript-r6.31.fasta ${references}/dm6/122_dmel_TE_Senti2015.fa) \
-i ${references}/dm6/indices/dm6_122TE_salmon_v1.1.0 -k 31


### common analyses for all RNA sequencing libraries --- START ---

### step 1: trim adapters, filter poor-quality reads, and take only the reads that have mates
### Libraries from this study and Senti 2015 used the same adapters.
mkdir -p ${analysis}/${lib_pA}/dm6_122TE_salmon_v1.1.0
### fastx_clipper from fastx_toolkit
fastx_clipper -l 35 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-i <(zcat ${raw_data}/${lib_pA}_R1.fastq.gz | fastq_quality_filter -q 33 -p 40) > ${analysis}/${lib_pA}/${lib_pA}_R1.trimmed.fastq
fastx_clipper -l 35 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGA \
-i <(zcat ${raw_data}/${lib_pA}_R2.fastq.gz | fastq_quality_filter -q 33 -p 40) > ${analysis}/${lib_pA}/${lib_pA}_R2.trimmed.fastq

### pair mates using fastq-pair and take only the reads that have mates
fastq_pair ${analysis}/${lib_pA}/${lib_pA}_R1.trimmed.fastq ${analysis}/${lib_pA}/${lib_pA}_R2.trimmed.fastq
for TYPE in R1 R2; do
gzip ${analysis}/${lib_pA}/${lib_pA}_${TYPE}.trimmed.fastq.paired.fq
rm ${analysis}/${lib_pA}/${lib_pA}_${TYPE}.trimmed.fastq.single.fq
rm ${analysis}/${lib_pA}/${lib_pA}_${TYPE}.trimmed.fastq
done

### step 2: mapping the paired-end reads to the dm6 transcriptome + 122 TEs using salmon
index_salmon="${references}/dm6/indices/dm6_122TE_salmon_v1.1.0"
salmon quant -i ${index_salmon} -l ISR --validateMappings --incompatPrior 0.0 --seqBias --gcBias \
-1 ${raw_data}/${lib_pA}_R1.trimmed.fastq.paired.fq.gz -2 ${raw_data}/${lib_pA}_R2.trimmed.fastq.paired.fq.gz \
-p 8 -o ${analysis}/${lib_pA}/dm6_122TE_salmon_v1.1.0/

### common analyses for all RNA sequencing libraries --- END ---

### ${scripts}/pA_RNAseq_salmon_DESeq2_TE.R was used to perform differential gene expression analysis. 
