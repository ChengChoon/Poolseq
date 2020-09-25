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

With the generated output file _rc: this file contains the major and minor alleles for every SNP in a concise format

We can calculate hetrozygosity for each pool seperatedly, or a pooled heterozygosity, Hp using the formula

Hp = 2*Sum(nMAJ)*Sum(nMIN)/[Sum(nMAJ) + Sum(nMIN)]^2, where Sum(nMAJ) and Sum(nMIN) are breedpool-specific sums of nMAJ and, respectively, nMIN counted at all SNPs in the window (if using sliding windows approach.

Please refer to python script, 3_pooled_het.py to do this. Where it based on python3 version. ref.bed (scaffolds-length) is required.

To prevent windows containing very few SNPs from adding spurious fixation signals, we ommited windows where only 1-10 SNPs had been detected.

If needed, we can Z-transform the Hp values to Zhp= (Hp-mean.Hp)/std_dev.Hp
 
## Calculate Fst for every SNP using fst-sliding.pl (also the parameter meaning) from popoolation2 

perl fst-sliding.pl --input p1_p2.sync --output p1_p2.fst --suppress-noninformative --min-count 6 --min-coverage 50 --max-coverage 200 --window-size 1 --step-size 1 --pool-size X:X

## Calculate Fst values using a sliding window approach

perl fst-sliding.pl --input p1_p2.sync --output p1_p2_w500.fst --min-count 6 --min-coverage 50 --max-coverage 200 --window-size 500 --step-size 250 --pool-size X:X

## plot the figure using Rscipt = plotscript.R or to rearrange the scaffolds into chromosome or LG groups. Please use forward.sh, reverse.sh, script.sh.

## reference from goat genome, https://doi.org/10.1371/journal.pgen.1008536

(i) Genome Analysis Toolkit UnifiedGenotyper version 3.7 with setting: -glm SNP, -stand_call_conf 20, -out_mode EMIT_VARIANTS_ONLY and ploidy 16/20/24, filtered using the generic hard-filtering recommendations available from https://gatkforums.broadinstitute.org/gatk/discussion/6925/understanding-and-adapting-thegeneric-hard-filtering-recommendations

(ii) SAMtools mpileup with the settings -q 15, -Q 20, -C 50 and -B, mpileup files were streamed to the PoPoolation2 used the scripts mpileup2sync.jar with settings --fastqtype sanger --min-qual 20, snp-frequency-diff.pl with the settings --min-coverage 15 --max-coverage 50 --min-count 3





