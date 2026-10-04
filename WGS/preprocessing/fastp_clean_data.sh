#!/bin/bash
# fastp_clean_data.sh
# Filter paired-end WGS raw data (.fq.gz) provided by the sequencing company.
#
# The filtering parameters (-g -q 5 -u 50 -n 15 -l 150 ...) were specified by
# the sequencing vendor and should be kept consistent across all samples.
#
# Usage:
#   1. Edit <sample_id> and <project_id> below.
#   2. chmod a+x fastp_clean_data.sh
#   3. nohup ./fastp_clean_data.sh > clean_data.log 2>&1 &
#
# Output: two cleaned .fq.gz files plus an HTML/JSON QC report per sample.

for name in CAP001                 # sample to process — edit per run
do
    project_id="<project_id>"      # sequencing project directory name
    cleaned_dir="/path/to/project/Cleaned_Data"

    input1="/path/to/project/${project_id}/01.RawData/${name}/${name}_1.fq.gz"
    input2="/path/to/project/${project_id}/01.RawData/${name}/${name}_2.fq.gz"
    output1="${cleaned_dir}/${name}_1.fq.gz"
    output2="${cleaned_dir}/${name}_2.fq.gz"
    html_rep="${cleaned_dir}/${name}.html"
    json_rep="${cleaned_dir}/${name}.json"

    fastp -i ${input1} -I ${input2} \
          -g -q 5 -u 50 -n 15 -l 150 \
          --overlap_diff_limit 1 --overlap_diff_percent_limit 10 \
          -o ${output1} -O ${output2} -h ${html_rep} -j ${json_rep}

    echo "${name} cleaned complete"
done
