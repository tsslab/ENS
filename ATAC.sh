#!/bin/sh

module load bedtools
module load bowtie/1.0.0
module load picard-tools/1.83
module load ucsctools
module load samtools
module load perl
module load macs2/2.0.10
module load bamtools
module load deeptools
module load python
module load basespace

CHROM=/t1-data/user/iling/genomes/galGal5.chrom.sizes
GENOME=/t1-data/user/iling/genomes/galGal5

AN=R1_Tpos

#To download samples from Basespace
samples2files -K {CLIENT_KEY} -S {CLIENT_SECRET} -A {CLIENT_ACCESS} -s {SAMPLE_ID}

gunzip *.gz

#Concatenate files
cat $AN*_R1_*.fastq > $AN\.R1.fastq
cat $AN*_R2_*.fastq > $AN\.R2.fastq

FA1=$AN\.R1.fastq
FA2=$AN\.R2.fastq

#Aligning sequence. the -X2000 ensures fragments up to 2kb were allowed to align and -m1 for unique reads only
bowtie -S -p 8 -m 1 -X 2000 $GENOME -1 $FA1 -2 $FA2  --chunkmb 500 $AN\.sam

#Get "smoothened" bigwig (using script from Jim Hughes). This will generate a smoothened BW file to be used on genome browser
~/scripts/sam2bwPE_gg5.pl -sam $AN\.sam -build galGal5 -name smo_$AN

#For downstream analysis and peak call
#Convert SAM to BAM file
samtools view -bS $AN\.sam > $AN\.bam

#Retain only aligned pairs based on BAM flag and remove chrM within the BAM file
samtools view -h $AN\.bam | awk '{ if ($3 != "chrM"){print $0; }}' |\
samtools view -Sb -b -F 4 - > $AN\.PAIRED_NOMT.bam

#sort BAM files by position
samtools sort $AN\.PAIRED_NOMT.bam -o $AN\.PAIRED_NOMT_nsort.bam

#remove duplicates
MarkDuplicates INPUT=$AN\.PAIRED_NOMT_nsort.bam OUTPUT=$AN\.NODUPS_NOMT.bam METRICS_FILE=$AN\.dup.txt \
REMOVE_DUPLICATES=True ASSUME_SORTED=True MAX_FILE_HANDLES_FOR_READ_ENDS_MAP= 1000

#Then sort by chromosome
samtools sort -n $AN\.NODUPS_NOMT.bam > $AN\.NODUPS_NOMT.nsort.bam

#Genome browser visualization
genomeCoverageBed -bg -split -ibam $AN\.NODUPS_NOMT.nsort.bam -g $CHROM > $AN\.bg
bedGraphToBigWig $AN\.bg $CHROM $AN\.bw

#Get fragment length
samtools view $AN\.NODUPS_NOMT.nsort.bam | awk '$9>0' | \
cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > $AN\.frag.txt

#Get BED file for peak call
#Then remove alignments that cause problems and then convert the BAM to BED file using -bedpe for
#paired end reads
samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 | \
bedtools bamtobed -bedpe -i $AN\.NODUPS_NOMT.nsort.bam > $AN\.NODUPS_NOMT.nsort.bed

#Pile up fragments and sort
cat $AN\.NODUPS_NOMT.nsort.bed | cut -f1,2,6,7,8,9,10 | \
sort -k1,1 -k2,2n > $AN\.NODUPS_NOMT_cut.bed

#Extending reads for MACS peak call 
bedtools slop -i $AN\.NODUPS_NOMT_cut.bed -g $CHROM -b 75 -s > $AN\.NODUPS_NOMT_cut_75shift.bed
