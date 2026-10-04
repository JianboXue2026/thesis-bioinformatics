#!/bin/bash
# samtools_fixmate.sh
# Fill in mate coordinates, ISIZE and mate-related flags for name-sorted
# alignments. This is a prerequisite for PCR duplicate marking with
# `samtools markdup`.
#
# Two variants are provided:
#   A) input BAM is name-sorted (produced by bwa_alignment.sh) — direct call
#   B) input BAM is coordinate-sorted — re-sort by name first (pipe)

core_processing=32

for name in CAP001   # sample to process — edit per run
do
    input_hg38="/path/to/project/BMA/HG38/${name}-hg38.bam"
    output_hg38="/path/to/project/BMA/HG38/${name}-hg38.sort.fixmate.bam"
    input_CHM13="/path/to/project/BMA/CHM13/${name}-CHM13.bam"
    output_CHM13="/path/to/project/BMA/CHM13/${name}-CHM13.sort.fixmate.bam"
    input_YAO="/path/to/project/BMA/T2T-YAO/${name}-YAO.bam"
    output_YAO="/path/to/project/BMA/T2T-YAO/${name}-YAO.sort.fixmate.bam"

    # --- A) input already name-sorted --------------------------------------
    samtools fixmate -@ ${core_processing} -m ${input_hg38} ${output_hg38}
    samtools fixmate -@ ${core_processing} -m ${input_CHM13} ${output_CHM13}
    samtools fixmate -@ ${core_processing} -m ${input_YAO} ${output_YAO}

    # --- B) input coordinate-sorted: re-sort by name, then fixmate ---------
    # samtools sort -@ ${core_processing} -T /path/to/temp/${name}-hg38.temp.fixmate -n -m 3G ${input_hg38} | \
    #     samtools fixmate -@ ${core_processing} -m - - > ${output_hg38}
    # (repeat for CHM13 / YAO accordingly)

    echo "${name} fixmate complete"
done
