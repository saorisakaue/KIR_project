#!/bin/bash

## HaplotypeCaller with ploidy per sample and per gene
java  -Xmx20g -Xms20g -jar /path/to/GATK3.X/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-nct 2 \
-R ./REF/KIR_seq_ref.fasta \
-I ${BAMPATH}/final/${ID}.mapped.rg.md.bam \
-o ${GVCFPATH}/${ID}.${GENE}.g.vcf.gz \
-ploidy ${PLOIDY} \
-L ./REF/${GENE}.intervals \
--emitRefConfidence GVCF
