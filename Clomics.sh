#!/bin/bash
# +
## SPECIFY:

## path to PHM
export clomics_dir=/software/clomics-1.0/bin/clomics;

## path to PHM helper scripts:
path_to_help=/data/pushkare/computational_paper/GITHUB/00.CM_mapping/cm_1.0-1_x86_64/usr/local/bin;

## input dataset:
dataset="test_data";

## path to gzipped BED count matrices
path_to_input=/data/pushkare/computational_paper/GITHUB/${dataset}/gzipped_bed;

## chromosomes to use
chromosomes="chr21,chr22";

## output path:
output_path=/data/pushkare/computational_paper/GITHUB/00.CM_mapping/Clomics;


## Run PHM on test_data for chr21 and chr22
sh ${path_to_help}/clomics/run_Clomics.sh \
    -i ${path_to_input}/H3K27ac_chr21_chr22.bed.gz,${path_to_input}/H3K4me1_chr21_chr22.bed.gz \
    -p ${clomics_dir} \
    -c ${chromosomes} \
    -s ${path_to_help} \
    -d ${dataset} \
    -o ${output_path}/${dataset}

# -


