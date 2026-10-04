#!/bin/bash
# bwa_alignment.sh
# Align paired-end clean reads to the three reference genomes with bwa mem,
# and pipe the output through `samtools sort -n` (name-sorted BAM), which is
# the sorting order required by the subsequent `samtools fixmate` step.
#
# The Read Group string (-R) is mandatory: GATK refuses BAM files without RG.
#
# Usage:
#   1. Edit the sample name and paths.
#   2. nohup ./bwa_alignment.sh > bwa-<sample>.log 2>&1 &
#
# Tip: inspect a resulting BAM (header, RG) with:
#   samtools view -h <sample>-YAO.bam | less -S

core_processing=16   # threads for bwa mem
core_converting=16   # threads for samtools sort

for name in CAP021   # sample to process — edit per run
do
    input_1="/path/to/project/Cleaned_Data/${name}_1.fq.gz"
    input_2="/path/to/project/Cleaned_Data/${name}_2.fq.gz"

    output_hg38="/path/to/project/BMA/HG38/${name}-hg38.bam"
    temp_dir_hg38="/path/to/temp/${name}-hg38.temp.bwa"
    output_CHM13="/path/to/project/BMA/CHM13/${name}-CHM13.bam"
    temp_dir_CHM13="/path/to/temp/${name}-CHM13.temp.bwa"
    output_YAO="/path/to/project/BMA/T2T-YAO/${name}-YAO.bam"
    temp_dir_YAO="/path/to/temp/${name}-YAO.temp.bwa"

    RG="@RG\tID:${name}\tSM:${name}\tLB:WES\tPL:ILLUMINA"

    bwa mem -t ${core_processing} -R ${RG} /path/to/ref_seq/hg38/hg38.fa ${input_1} ${input_2} | \
        samtools sort -n -m 3G -@ ${core_converting} -T ${temp_dir_hg38} -o ${output_hg38} -  &

    bwa mem -t ${core_processing} -R ${RG} /path/to/ref_seq/CHM13/T2T-CHM13v2.0.fna ${input_1} ${input_2} | \
        samtools sort -n -m 3G -@ ${core_converting} -T ${temp_dir_CHM13} -o ${output_CHM13} -  &

    bwa mem -t ${core_processing} -R ${RG} /path/to/ref_seq/YAO/T2T-Yao-hp.v1.1.fasta ${input_1} ${input_2} | \
        samtools sort -n -m 3G -@ ${core_converting} -T ${temp_dir_YAO} -o ${output_YAO} -  &
done
