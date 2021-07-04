#!/bin/bash

VARIANT=`ls ${GVCFPATH}/ | grep -v tbi | grep ${GENE} | awk -v p=${GVCFPATH} '{printf("--variant %s/%s ",p,$1)}'`

java -Xmx40g -Xms40g -jar /path/to/GATK3.X/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R ./REF/KIR_seq_ref.fasta \
-allSites \
-o ${VCFPATH}/${GENE}.JointGenotype.vcf.gz \
${VARIANT}
