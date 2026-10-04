#!/bin/bash
# bwa_index.sh
# Build BWA indexes for the three reference genomes used in this study:
#   hg38 (GRCh38), T2T-CHM13 and T2T-YAO.
#
# Each genome directory ends up with a set of index files prefixed with the
# genome file name (.amb / .ann / .bwt / .pac / .sa).

# hg38
cd /path/to/ref_seq/hg38/
bwa index -a bwtsw hg38.fa

# CHM13
cd /path/to/ref_seq/CHM13/
bwa index -a bwtsw T2T-CHM13v2.0.fna

# T2T-YAO
cd /path/to/ref_seq/YAO/
bwa index -a bwtsw T2T-Yao-hp.v1.1.fasta
