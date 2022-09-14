#!/bin/bash

SAMPLE=$1

REFERENCE_FASTA="/Users/mgjones/projects/ecDNA_tracing/REF/ref_${SAMPLE}_long.fa"

FASTQ_HOME="/Users/mgjones/projects/ecDNA_tracing/SEQ/pilot3/nextseq/"
SAMPLE_R1="${FASTQ_HOME}/${SAMPLE}"*R1*.fastq.gz
SAMPLE_R2="${FASTQ_HOME}/${SAMPLE}"*R2*.fastq.gz

MERGED_R1="${FASTQ_HOME}/${SAMPLE}_R1_merged.fastq.gz"
MERGED_R2="${FASTQ_HOME}/${SAMPLE}_R2_merged.fastq.gz"

INDEX="/Users/mgjones/projects/ecDNA_tracing/bowtie_index/ptracer_${SAMPLE}_long"
OUTPUT_DIRECTORY="/Users/mgjones/projects/ecDNA_tracing/mapped/pilot3/nextseq/"
SAMPLE_MAPPED="${OUTPUT_DIRECTORY}/${SAMPLE}_edited.sam"

PICARD="/Users/mgjones/software/picard.jar"

# Specify number of threads
THREADS=10

# Specify min quality
MIN_QUALITY=25

# set up directories
mkdir -p $OUTPUT_DIRECTORY

# concatenate R1 and R2
echo ">> Concatenating reads..."
# cat $SAMPLE_R1 > $MERGED_R1
# cat $SAMPLE_R2 > $MERGED_R2

echo ">> Masking low-quality bases..."
# python /Users/mgjones/projects/ecDNA_tracing/scripts/quality_control_reads.py \
#     $MERGED_R1 $MERGED_R2 --min_quality $MIN_QUALITY

echo ">> Compressing new masked fastqs..."
# MASKED_R1="${MERGED_R1%.fastq.gz}".masked.fastq
# MASKED_R2="${MERGED_R2%.fastq.gz}".masked.fastq
# pigz -p 10 $MASKED_R1
# pigz -p 10 $MASKED_R2

# create Bowtie indices
if [ ! -f "$INDEX.1.bt2" ]; then
    echo ">> Building Reference..."
    bowtie2-build $REFERENCE_FASTA $INDEX
fi

# map reads
echo ">> Aligning reads..."
MASKED_R1="${MERGED_R1%.fastq.gz}".masked.fastq.gz
MASKED_R2="${MERGED_R2%.fastq.gz}".masked.fastq.gz
bowtie2 -x $INDEX -1 $MASKED_R1 -2 $MASKED_R2 -S $SAMPLE_MAPPED \
    -p $THREADS --very-sensitive -X 700 \
    --n-ceil 3 \
    --trim3 10 --trim5 10 --rg-id "ID:${SAMPLE}\tSM:${SAMPLE}\tLB:None\tPL:Illumina"

# sort reads
echo ">> Sorting alignments..."
SAMPLE_BAM_SORTED="${SAMPLE_MAPPED%.sam}".sorted.bam
samtools sort -@ $THREADS $SAMPLE_MAPPED -o $SAMPLE_BAM_SORTED

# reindex
samtools index -@ $THREADS $SAMPLE_BAM_SORTED

# filter out PCR duplicates
# echo ">> Removing duplicates..."
# SAMPLE_RMDUP="${SAMPLE_MAPPED%.sam}".rmdup.bam
# java -jar $PICARD MarkDuplicates INPUT=$SAMPLE_BAM_SORTED OUTPUT=$SAMPLE_RMDUP \
#     METRICS_FILE="${OUTPUT_DIRECTORY}/${SAMPlE}".PicardMetrics.txt \
#     VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true > "${OUTPUT_DIRECTORY}/${SAMPLE}".picard.log

# filter reads by length and MAPQ
echo ">> Filtering alignments..."
SAMPLE_BAM_FILTERED="${SAMPLE_MAPPED%.sam}".filtered.bam
samtools view -h $SAMPLE_BAM_SORTED | \
    awk 'length($10) > 40 || $1 ~ /^@/' | \
    samtools view -bSq 30 > $SAMPLE_BAM_FILTERED

# reindex
samtools index -@ $THREADS $SAMPLE_BAM_FILTERED
