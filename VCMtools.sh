#!/bin/bash
# +
## SPECIFY:

## path to VCMtools helper scripts:
path_to_help=./vcmtools;

## input dataset:
dataset="test_data";

## path to gzipped BED count matrices
path_to_input=./${dataset}

## chromosomes to use
chromosomes="chr21,chr22";
pvalue=0.001;
n_cores=1;

## output path:
output_path=./VCMtools;

## Run VCMtools on test_data for chr21 and chr22
sh ${path_to_help}/vcmtools/run_VCMtools.sh \
    -s ${path_to_help} \
    -d ${dataset} \
    -f ${path_to_input}/H3K27ac_chr21_chr22.bed,${path_to_input}/H3K4me1_chr21_chr22.bed \
    -c ${chromosomes} \
    -p ${pvalue} \
    -n ${n_cores} \
    -o ${output_path}/${dataset}
# -


