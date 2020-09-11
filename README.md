# Poolseq

## Same procedures as GWAS for the raw data processing steps:

#concatenate several reads fastq files together

#trim ends or adapters using trimmomatic-0.39.jar

#check the output of trimmed files with fastQC:

#mapping using bwa mem

#convert file format for size reduction sam to bam

#certain program required sorted bam file for the following steps

#create index for Picard to run

#run Bedtools to check the overal genomic coverage of the reads

#run Picard to mark & remove duplicated reads

## filter reads with quality lower than 20 & flag 1804 (check manual of samtools)

samtools view -@ 16 -q 20 -F 1804 -b R.sorted.MD.bam > R.sorted.MD.Q20.bam 

## combine several groups of pool into a file

samtools mpileup -B -f ref.fa pool1.bam pool2.bam > p1.p2.mpileup

## synchronizing the mpileup file using java script mpileup2sync.jar from popoolation2

java -Xmx32G -jar mpileup2sync.jar --input p1.p2.mpileup --output p1.p2.java.sync --fastq-type sanger --min-qual 20 --threads 8

## Calculate allele frequency differences using snp-frequency-diff.pl (also the parameter meaning) from popoolation2 

perl snp-frequency-diff.pl --input p1_p2.java.sync --output-prefix p1_p2 --min-count 6 --min-coverage 50 --max-coverage 200

## Calculate Fst for every SNP using fst-sliding.pl (also the parameter meaning) from popoolation2 

perl fst-sliding.pl --input p1_p2.sync --output p1_p2.fst --suppress-noninformative --min-count 6 --min-coverage 50 --max-coverage 200 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size X:X

## Calculate Fst values using a sliding window approach

perl fst-sliding.pl --input p1_p2.sync --output p1_p2_w500.fst --min-count 6 --min-coverage 50 --max-coverage 200 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size X:X

## plot the figure using Rscipt = plotscript.R
