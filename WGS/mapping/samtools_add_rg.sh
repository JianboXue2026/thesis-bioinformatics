#!/bin/bash
# samtools_add_rg.sh
# Rescue script: re-add the Read Group (RG) to BAM files that were aligned
# without `-R`, which would otherwise be rejected by GATK downstream.

for name in CAP001   # sample to process — edit per run
do
    RG="@RG\tID:${name}\tSM:${name}\tLB:WES\tPL:ILLUMINA"

    input_hg38="/path/to/project/BMA/HG38/${name}.bam"
    output_hg38="/path/to/project/BMA/HG38/${name}-hg38.bam"
    input_CHM13="/path/to/project/BMA/CHM13/${name}.bam"
    output_CHM13="/path/to/project/BMA/CHM13/${name}-CHM13.bam"
    input_YAO="/path/to/project/BMA/T2T-YAO/${name}.bam"
    output_YAO="/path/to/project/BMA/T2T-YAO/${name}-YAO.bam"

    samtools addreplacerg -r ${RG} -@ 4 -o ${output_hg38}  ${input_hg38}
    samtools addreplacerg -r ${RG} -@ 4 -o ${output_CHM13} ${input_CHM13}
    samtools addreplacerg -r ${RG} -@ 4 -o ${output_YAO}   ${input_YAO}

    echo "${name} RG adding complete"
done
