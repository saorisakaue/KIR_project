#!/bin/bash

# align to KIR reference with bwa mem
bwa mem -t 2 \
./REF/KIR_seq_ref \
-R "@RG\tID:${ID}\tSM:${ID}" \
${SAMPLEPATH}_R1.fastq.gz \
${SAMPLEPATH}_R2.fastq.gz \
| samtools view -@ 2 -b - | samtools sort -@ 2 - -o ${BAMPATH}/raw/${ID}.mapped.bam

## Add ReadGroups
java -Xmx40g -Xms40g -jar /path/to/picard/picard.jar AddOrReplaceReadGroups \
I=${BAMPATH}/raw/${ID}.mapped.bam \
O=${BAMPATH}/interval/${ID}.mapped.rg.bam \
RGLB=${ID} RGPL=ILLUMINA RGPU=${ID} RGSM=${ID} RGID=${ID} \
VALIDATION_STRINGENCY=LENIENT

## Mark Duplicates
java -Xmx40g -Xms40g -jar /path/to/picard/picard.jar MarkDuplicates \
I=${BAMPATH}/interval/${ID}.mapped.rg.bam \
O=${BAMPATH}/final/${ID}.mapped.rg.md.bam \
ASSUME_SORTED=false \
REMOVE_DUPLICATES=false \
CREATE_INDEX=True \
VALIDATION_STRINGENCY=LENIENT \
M=${BAMPATH}/interval/${ID}.mapped.rg.md.metrics

## Diagnose targets
java -Xmx40g -Xms40g -jar /path/to/GATK3.X/GenomeAnalysisTK.jar \
-T DiagnoseTargets \
-R ./REF/KIR_seq_ref.fasta \
-o ${COVERAGE}/DiagnoseTargets.${ID}.coverage.vcf \
-I ${BAMPATH}/final/${ID}.mapped.rg.md.bam \
-L ./REF/KIR_seq_ref.intervals

