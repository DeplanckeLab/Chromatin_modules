#!/bin/bash
## Specify output folder
output_path="./clomics_test";

## Specify paths to Clomics executable and helper scripts
path_to_clomics="./clomics-2.0/bin/clomics";
path_to_main="./clomics/main";

## define Clomics-specific parameters
bg_threshold="3";  ## default bg_threshold="3"
n_peaks="200";  ## default n_peaks="200"

## Specify the number of processs to use
## if parallelizing CRD to CM conversion step
n_cores=1;  ## default n_cores=1

## Specify chromosomes to use for CM mapping
chromosomes="22";  ## dash defines a range, default is "1-22" (all chromosomes from chr1 to chr22);
## Alternatively, "chromosomes" parameter can be defined as
## "1,5,10" (without spaces) when using, in this example, three chromosomes;
## or "22", when using a single chromosome

## Specify cell type/dataset name
dataset_name="test_data_chr22";

## Specify paths to count matrices for ChIP-seq/ATAC-seq data
input_path="./test_data/gzipped_bed";
k4me1_matrix="H3K4me1:${input_path}/H3K4me1_chr22.bed.gz";
k27ac_matrix="H3K27ac:${input_path}/H3K27ac_chr22.bed.gz";

## Run Clomics pipeline with the variables defined above

Rscript ${path_to_main}/main_clomics_script.R \
    -path_to_clomics ${path_to_clomics} \
    -path_to_src "${path_to_main}/src" \
    -bg_threshold ${bg_threshold} \
    -n_peaks ${n_peaks} \
    -n_cores ${n_cores} \
    -chromosomes ${chromosomes} \
    -dataset_name ${dataset_name} \
    -matrix_files "${k4me1_matrix},${k27ac_matrix}" \
    -output_folder ${output_path}
