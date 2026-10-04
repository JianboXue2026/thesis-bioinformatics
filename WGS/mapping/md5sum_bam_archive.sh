#!/bin/bash
# md5sum_bam_archive.sh
# Generate md5 checksums for the pre-variant-calling BAM archive (BAM + BAM.BAI
# per sample and reference), so the archive can be verified after transfer to
# backup storage.
#
# Usage:
#   1. Edit start_num / end_num to the sample range.
#   2. Run inside the archive root (expects ./HG38/ ./CHM13/ ./T2T-YAO/).
#   3. Verify on the backup copy with:
#        nohup md5sum -c BMA_CAPxxx-CAPxxx.md5.txt >> BMA_CAPxxx-CAPxxx_md5-check.txt &

cd /path/to/project/BMA

start_num=1      # first sample number
end_num=20       # last sample number
start_num_formated=$(printf "%03d" ${start_num})
end_num_formated=$(printf "%03d" ${end_num})

prefix="CAP"
suffix_1=".bam"
suffix_2=".sort.fixmate.samtools_dedup.bam"
suffix_3=".sort.fixmate.samtools_dedup.bam.bai"

number_list="$(seq -f "${prefix}%03g" ${start_num} ${end_num})"

for ref_flag in "hg38" "CHM13" "T2T-YAO"
do
  # NOTE: the spaces around [ ] and = are required in bash tests
  if [ ${ref_flag} = "hg38" ]; then
    dir="./HG38/"
    midfix="-hg38"
  elif [ ${ref_flag} = "CHM13" ]; then
    dir="./CHM13/"
    midfix="-CHM13"
  else
    dir="./T2T-YAO/"
    midfix="-YAO"
  fi

  for name in ${number_list}
  do
    file_name_1="${dir}${name}${midfix}${suffix_1}"
    file_name_2="${dir}${name}${midfix}${suffix_2}"
    file_name_3="${dir}${name}${midfix}${suffix_3}"
    md5sum ${file_name_1} ${file_name_2} ${file_name_3} \
        >> "./BMA_CAP${start_num_formated}-CAP${end_num_formated}.md5.txt"
  done
done
