#!/bin/bash
# +
## SPECIFY:

## path to Clomics
export clomics_dir=/software/clomics-1.0/bin/clomics;

## path to Clomics helper scripts:
path_to_help=./;

## input dataset:
dataset="test_data";

## path to gzipped BED count matrices
path_to_input=./${dataset}/gzipped_bed;

## chromosomes to use
chromosomes="chr21,chr22";

## output path:
output_path=./Clomics;


## Run Clomics on test_data for chr21 and chr22
sh ${path_to_help}/clomics/run_Clomics.sh \
    -i ${path_to_input}/H3K27ac_chr21_chr22.bed.gz,${path_to_input}/H3K4me1_chr21_chr22.bed.gz \
    -p ${clomics_dir} \
    -c ${chromosomes} \
    -s ${path_to_help} \
    -d ${dataset} \
    -o ${output_path}/${dataset}

# -


