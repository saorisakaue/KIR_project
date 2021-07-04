#!/bin/bash

singularity exec \
 /path/to/deepvariant_0.10.0.simg /opt/deepvariant/bin/run_deepvariant \
 --model_type=WGS \
 --ref /path/to/REF/KIR_seq_ref.fasta \
 --reads ${BAMPATH}/final/${ID}.mapped.rg.md.bam \
 --output_vcf=${DeepVariantPATH}/${ID}.dv.vcf.gz \
 --output_gvcf=${DeepVariantPATH}/${ID}.dv.g.vcf.gz \
 --num_shards=2