#!/bin/bash
# +
path_to_scripts=/data/pushkare/computational_paper/GITHUB/00.CM_mapping/Clomics/main;
path_to_gz_bed=/data/pushkare/computational_paper/GITHUB/test_data/gzipped_bed;
output_path=/data/pushkare/computational_paper/GITHUB/00.CM_mapping/Clomics;

sh ${path_to_scripts}/Clomics.sh \
    -i ${path_to_gz_bed}/H3K27ac_chr21_chr22.bed.gz,${path_to_gz_bed}/H3K4me1_chr21_chr22.bed.gz \
    -s ${path_to_scripts} \
    -d "test_data" \
    -m "all_marks"\
    -o ${output_path}
# -


