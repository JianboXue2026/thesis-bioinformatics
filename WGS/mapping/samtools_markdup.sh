#!/bin/bash
# samtools_markdup.sh
# Remove PCR duplicates from fixmate-processed BAM files:
# coordinate-sort the fixmate output, then mark duplicates with samtools markdup.

core_processing=16

for name in CAP001   # sample to process — edit per run
do
    input_hg38="/path/to/project/BMA/HG38/${name}-hg38.sort.fixmate.bam"
    temp_dir_hg38="/path/to/temp/${name}-hg38.temp.dedup"
    output_hg38="/path/to/project/BMA/HG38/${name}-hg38.sort.fixmate.samtools_dedup.bam"

    input_CHM13="/path/to/project/BMA/CHM13/${name}-CHM13.sort.fixmate.bam"
    temp_dir_CHM13="/path/to/temp/${name}-CHM13.temp.dedup"
    output_CHM13="/path/to/project/BMA/CHM13/${name}-CHM13.sort.fixmate.samtools_dedup.bam"

    input_YAO="/path/to/project/BMA/T2T-YAO/${name}-YAO.sort.fixmate.bam"
    temp_dir_YAO="/path/to/temp/${name}-YAO.temp.dedup"
    output_YAO="/path/to/project/BMA/T2T-YAO/${name}-YAO.sort.fixmate.samtools_dedup.bam"

    samtools sort -@ ${core_processing} -T ${temp_dir_hg38} -m 3G ${input_hg38} | \
        samtools markdup -@ ${core_processing} - ${output_hg38}

    samtools sort -@ ${core_processing} -T ${temp_dir_CHM13} -m 3G ${input_CHM13} | \
        samtools markdup -@ ${core_processing} - ${output_CHM13}

    samtools sort -@ ${core_processing} -T ${temp_dir_YAO} -m 3G ${input_YAO} | \
        samtools markdup -@ ${core_processing} - ${output_YAO}

    echo "${name} samtools deduplication complete"
done
