#!/bin/bash

##### USAGE #####
# bash reference_based_pipeline.sh <R1.fastq.gz> <R2.fastq.gz> <reference.fa> <sample_name>

set -e

R1=$1
R2=$2
REF=$3
SAMPLE=$4

THREADS=16

echo "Starting pipeline for ${SAMPLE}"

########################################
# Step 1: Adapter Trimming
########################################
echo "Step 1: Adapter Trimming"

fastq-mcf -x 10 -q 20 -l 30 -S -k 0 \
-o ${SAMPLE}.AT.R1.fastq.gz \
-o ${SAMPLE}.AT.R2.fastq.gz \
adapters/adapters.fa $R1 $R2

########################################
# Step 2: Host (Human) Contamination Removal
########################################
echo "Step 2: Host Removal"

bwa mem -t $THREADS reference/hg19.fa \
${SAMPLE}.AT.R1.fastq.gz ${SAMPLE}.AT.R2.fastq.gz > ${SAMPLE}.con.sam

samtools view -bS -@ 8 ${SAMPLE}.con.sam > ${SAMPLE}.con.bam
samtools sort -@ 8 ${SAMPLE}.con.bam -o ${SAMPLE}.con.sorted.bam
samtools index ${SAMPLE}.con.sorted.bam

samtools flagstat ${SAMPLE}.con.sorted.bam > ${SAMPLE}.con.alnstats.txt

# Extract unmapped reads
samtools view -h -f 4 ${SAMPLE}.con.sorted.bam > ${SAMPLE}.unmapped.sam
samtools view -bS ${SAMPLE}.unmapped.sam > ${SAMPLE}.unmapped.bam
samtools sort -n ${SAMPLE}.unmapped.bam -o ${SAMPLE}.unmapped.sorted.bam

bedtools bamtofastq \
-i ${SAMPLE}.unmapped.sorted.bam \
-fq ${SAMPLE}.UA.R1.fastq.gz \
-fq2 ${SAMPLE}.UA.R2.fastq.gz

########################################
# Step 3: Reference Alignment
########################################
echo "Step 3: Reference Alignment"

bwa mem -t $THREADS $REF \
${SAMPLE}.UA.R1.fastq.gz ${SAMPLE}.UA.R2.fastq.gz > ${SAMPLE}.final.sam

samtools view -bS -@ 8 ${SAMPLE}.final.sam > ${SAMPLE}.final.bam
samtools sort -@ 8 ${SAMPLE}.final.bam -o ${SAMPLE}.final.sorted.bam
samtools index ${SAMPLE}.final.sorted.bam

########################################
# Step 4: Add Read Groups
########################################
echo "Step 4: Adding Read Groups"

picard AddOrReplaceReadGroups \
I=${SAMPLE}.final.sorted.bam \
O=${SAMPLE}.final.sorted.withtags.bam \
RGID=${SAMPLE} \
RGLB=LIB \
RGPL=illumina \
RGSM=${SAMPLE} \
RGPU=unit1 \
CREATE_INDEX=true

samtools flagstat ${SAMPLE}.final.sorted.withtags.bam > ${SAMPLE}.final.alnstats.txt

########################################
# Step 5: Consensus Generation
########################################
echo "Step 5: Consensus Generation"

samtools mpileup -uf $REF ${SAMPLE}.final.sorted.withtags.bam | \
bcftools call -c | \
vcfutils.pl vcf2fq > ${SAMPLE}.fastq

seqtk seq -aQ64 -q20 -n N ${SAMPLE}.fastq > ${SAMPLE}_consensus.fasta

########################################
# Step 6: Coverage Statistics
########################################
echo "Step 6: Coverage Analysis"

samtools faidx $REF
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${REF}.fai > reference.bed

bamdst -p reference.bed -o coverage_${SAMPLE} ${SAMPLE}.final.sorted.withtags.bam

########################################
# Step 7: Variant Calling
########################################
echo "Step 7: Variant Calling"

gatk HaplotypeCaller \
-R $REF \
-I ${SAMPLE}.final.sorted.withtags.bam \
-O ${SAMPLE}_variants.vcf

echo "Pipeline completed for ${SAMPLE}"
