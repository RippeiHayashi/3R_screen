### the code was used to analyse the small RNA sequencing data for ovary samples from Drosophila melanogaster.

### software used
### fastx_toolkit
http://hannonlab.cshl.edu/fastx_toolkit/
### bowtie-1.2.3-linux-x86_64
https://github.com/BenLangmead/bowtie/releases/tag/v1.2.3
### samtools/1.19
https://github.com/samtools/
### bedtools/2.28.0
https://github.com/arq5x/bedtools2/releases/tag/v2.28.0
### kent-ucsc tools
https://github.com/ucscGenomeBrowser/kent
### R/4.0.0
https://cran.r-project.org/src/base/R-4/R-4.0.0.tar.gz
### weblogo 3.7.8
https://github.com/WebLogo/weblogo

### define the variables
scripts="<directory where the scripts are kept>"
references="<directory where the genome sequence file and other sequence files are kept>"
raw_data="<directory where the fastq files are kept>"
analysis="<directory where the output files are kept>"
lib_sRNA="<name of the small RNA seq library>"
### directory to save plots per library
mkdir -p ${analysis}/${lib_sRNA}/plots


### processing and analysing the small RNA sequencing libraries --- START ---

### below are the small RNA libraries processed in this script:
w1118_ovaries_sRNA_rep1
w1118_ovaries_sRNA_rep2
spnE_ZF_mutant_ovaries_sRNA_rep1
spnE_ZF_mutant_ovaries_sRNA_rep2
spnE_null_mutant_ovaries_sRNA_rep1
spnE_null_mutant_ovaries_sRNA_rep2


### preparing for the analyses, do this before processing individual libraries --- START ---
### step 0-1: generating a bowtie index for miscellaneous RNAs
### step 0-1-1: download miscRNA sequences from the FlyBase
mkdir -p ${references}/dm6/indices
wget --directory-prefix="${references}/dm6" http://ftp.flybase.net/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-miscRNA-r6.31.fasta.gz
wget --directory-prefix="${references}/dm6" http://ftp.flybase.net/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-tRNA-r6.31.fasta.gz
wget --directory-prefix="${references}/dm6" http://ftp.flybase.net/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-ncRNA-r6.31.fasta.gz

### step 0-1-2: making fasta including miRNA, rRNA, snRNA, snoRNA, and tRNA sequences
zcat ${references}/dm6/dmel-all-miscRNA-r6.31.fasta.gz ${references}/dm6/dmel-all-miRNA-r6.31.fasta.gz ${references}/dm6/dmel-all-tRNA-r6.31.fasta.gz |\
fasta_formatter - -t | tr -d ";" |\
awk '{split($5,a,"="); if($2~"miRNA" || $2~"rRNA" || $2~"snRNA" || $2~"snoRNA" || $2~"tRNA") print ">"a[2]"\n"$NF
}' > ${references}/dm6/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31.fasta

### step 0-1-3: making a bowtie index
bowtie-build ${references}/dm6/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31.fasta ${references}/dm6/indices/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31

### step 0-2: generating a bowtie index for endo siRNAs
zcat ${references}/dm6/dmel-all-ncRNA-r6.31.fasta.gz | fasta_formatter - -t | awk '{if($0 ~ "hpRNA") print ">"$1":hpRNA""\n"$NF}' > ${references}/dm6/dmel-hpRNA-r6.31.fasta
bowtie-build ${references}/dm6/dmel-hpRNA-r6.31.fasta ${references}/dm6/indices/dmel-hpRNA-r6.31

### step 0-3: make a bowtie index for the Drosophila melanogaster transposon sequences used in Senti G&D 2015, PMID: 26302790 (122 entries)
bowtie-build ${references}/TEs/122_dmel_TE_Senti2015.fa ${references}/TEs/indices/122_dmel_TE_Senti2015
fasta_formatter -i ${references}/TEs/122_dmel_TE_Senti2015.fa -t | awk '{print $1,length($2)}' > ${references}/TEs/122_dmel_TE_Senti2015.sizes

### step 0-4: generate pseudocounts of 1 at every coordinate of non-redundant TEs
mkdir -p ${references}/TEs/TE_pseudo.piRNAs/
cat ${references}/TEs/122_dmel_TE_Senti2015.sizes | while read TE length; do
if [[ ! -f ${references}/TEs/TE_pseudo.piRNAs/${TE}_pseudo.piRNAs.txt ]]; then
### generate pseudocounts of 1 at every coordinate
for ((i=0; i<=${length}-1; i++)); do
j=$((${i}+1))
printf "${TE} ${i} ${j} AAAAAAAAAAAAAAAAAAAAAAA@1 . +""\n""${TE} ${i} ${j} AAAAAAAAAAAAAAAAAAAAAAA@1 . -""\n" |\
tr ' ' '\t' >> ${references}/TEs/TE_pseudo.piRNAs/${TE}_pseudo.piRNAs.txt
done
fi

### common analyses for all small RNA libraries --- START ---

### step 1: trim the adapter sequence and collapse reads by sequence
### fastq_to_fasta, fasta_formatter, and fastx_clipper are from fastx_toolkit
fastq_to_fasta -Q33 -i <(zcat ${raw_data}/${lib_sRNA}.fastq.gz) |\
# -c: discard non-clipped sequences, -l 18: discard sequences shorter than 18nt
# collapse reads of the same sequence to one entry while retaining the number of reads
fastx_clipper -c -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -l 18 -i - | fasta_formatter - -t |\
awk '{if(length($NF)>25 && length($NF)<49) READS[substr($NF,5,length($NF)-8)]++} END {
for(var in READS) print ">"var"@"READS[var]"\n"var}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa

### ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa looks as follows:
>TCCAGATGAACGGTTAAGTGTCCAAAAAG@12
TCCAGATGAACGGTTAAGTGTCCAAAAAG
>TACTTGAAAGAATCAGGGGCCAACCAT@5
TACTTGAAAGAATCAGGGGCCAACCAT
...


### step 2: map to the Drosophila melanogaster miscellaneous RNA (misc RNA) allowing up to 1MM
### run bowtie to map reads to the miscRNA and take unmapped reads
index_misc="${references}/dm6/indices/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31"
bowtie -f -v 1 -a -S ${index_misc} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa \
--un ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-mapped.bed


### step 3: map to the Drosophila melanogaster endogenous siRNAs allowing up to 1MM
index_endo_siRNAs="${references}/dm6/indices/dmel-hpRNA-r6.31"
bowtie -f -v 1 --all --best --strata -S ${index_endo_siRNAs} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped_endo_siRNAs.bed


### step 4: mapping and characterisation of transposon small RNAs
### step 4-1: bowtie mapping to TEs, --all --best --strata option, allowing up to three mismatches
index_TE="${references}/TEs/indices/122_dmel_TE_Senti2015"
bowtie -f -v 3 --all --best --strata -S ${index_TE} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -b -S - | bedtools bamtobed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.bed

### step 4-2: plotting 5' and 3' ends of piRNAs mapping to TEs, and measure the ping-pong signature --- START ---
mkdir -p ${analysis}/${lib_sRNA}/End.plots.TEs

### step 4-2-1: extract top 120 TEs (sense and antisense) that produce abundant piRNAs, which gives ~ 80 top TEs
awk '{split($4,a,"@"); if($3-$2>22) TE[$1" "$6]+=a[2]} END {for(var in TE) print var,TE[var]}' ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.bed |\
sort -k3,3nr | head -n 120 > ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.top120

### run the code per TE --- START ---
awk '{print $1}' ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.top120 | sort | uniq | while read TE; do
if [[ ! -f ${analysis}/${lib_sRNA}/End.plots.TEs/${lib_sRNA}_${TE}_3MM_mappers.pingpong.plus_5end.minus_5end.table ]]; then
### count 5end and 3end of plus and minus mapping piRNA reads (>22nt)
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.bed ${references}/TEs/TE_pseudo.piRNAs/${TE}_pseudo.piRNAs.txt |\
awk -v TE=${TE} '{split($4,a,"@"); if($6=="+" && length(a[1])>22 && $1==TE) {PLUS5[$2]+=a[2]; PLUS3[$3-1]+=a[2]
} else if($6=="-" && length(a[1])>22 && $1==TE) {MINUS5[$3-1]+=a[2]; MINUS3[$2]+=a[2]
}} END {for(var in PLUS5) print var,PLUS5[var]-1,"plus_5end""\n"var,PLUS3[var]-1,"plus_3end""\n"var,MINUS5[var]-1,"minus_5end""\n"var,MINUS3[var]-1,"minus_3end"
}' > ${analysis}/${lib_sRNA}/End.plots.TEs/${lib_sRNA}_${TE}_3MM_mappers.counts

DIRECTORY="${analysis}/${lib_sRNA}/End.plots.TEs"
length=`(awk -v TE=${TE} '{if($1==TE) print $2}' ${references}/TEs/122_dmel_TE_Senti2015.sizes)`
### plotting 5' and 3' ends of sense and antisense mapping piRNAs
Rscript ${scripts}/plotting_TE-piRNAs.R DIRECTORY=${DIRECTORY} LIB=${lib_sRNA} TE=${TE} LENGTH=${length}
### measure pingpong signature
Rscript ${scripts}/pingpong_linkage.R DIRECTORY=${DIRECTORY} LIB=${lib_sRNA} TE=${TE}
fi
done
### run the code per TE --- END ---


### step 5: use weblogo to look at the nucleotide composition around the 3'end of TE antisense piRNAs --- START ---
mkdir -p ${analysis}/${lib_sRNA}/weblogo

TE_fasta="${references}/TEs/122_dmel_TE_Senti2015.fa"
cat ${references}/TEs/122_dmel_TE_Senti2015.sizes | while read TE length; do
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.bed | awk -v TE=${TE} '{if($1==TE && $3-$2>22) print}' |\
awk '{split($4,a,"@"); split(a[2],b,":"); for(i=1;i<=b[1];i++) {if($6=="-") print $1,$2-6,$2+5,$4,$5,$6}}' |\
awk -v LENGTH=${length} '{if($2>=0 && $3<=LENGTH) print}' | tr ' ' '\t' |\
bedtools getfasta -s -fi ${TE_fasta} -tab -bed - | awk '{print ">"$1"\n"toupper($2)}' | tr 'T' 'U' >> ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_all-TE-AS-piRNAs.3end_11nt_window.fasta
done

TYPE="all-TE-AS-piRNAs.3end_11nt_window"
weblogo -U probability -A rna -f ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TYPE}.fasta -F pdf -n 50 -c classic -o ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TYPE}.gt22_logo_prob.pdf
weblogo -U probability -A rna -f ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TYPE}.fasta -F logodata -n 50 -c classic -s large -t ${lib_sRNA}"_"${TE} -o ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TYPE}.gt22_logo_prob.txt

### common analyses for all small RNA libraries --- END ---
