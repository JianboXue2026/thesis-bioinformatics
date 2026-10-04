#!/bin/bash
# samtools_index_stats.sh
# Index the deduplicated BAM files and generate `samtools stats` reports.
# The resulting <sample>-<ref>-stat.txt files are summarised by
# preprocessing/extract_bamstat_info.py.

core_processing=8

for name in CAP001   # sample to process — edit per run
do
    bam_hg38="/path/to/project/BMA/HG38/${name}-hg38.sort.fixmate.samtools_dedup.bam"
    bam_CHM13="/path/to/project/BMA/CHM13/${name}-CHM13.sort.fixmate.samtools_dedup.bam"
    bam_YAO="/path/to/project/BMA/T2T-YAO/${name}-YAO.sort.fixmate.samtools_dedup.bam"
    stats_dir="/path/to/project/Samtools_stats"

    samtools index -@ ${core_processing} ${bam_hg38}  ${bam_hg38}.bai
    samtools index -@ ${core_processing} ${bam_CHM13} ${bam_CHM13}.bai
    samtools index -@ ${core_processing} ${bam_YAO}   ${bam_YAO}.bai

    samtools stats -@ ${core_processing} ${bam_hg38}  > "${stats_dir}/${name}-hg38-stat.txt"
    samtools stats -@ ${core_processing} ${bam_CHM13} > "${stats_dir}/${name}-CHM13-stat.txt"
    samtools stats -@ ${core_processing} ${bam_YAO}   > "${stats_dir}/${name}-YAO-stat.txt"

    echo "${name} stats complete"
done
