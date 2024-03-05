#!/bin/bash
Help()
{
   # Display Help
   echo "Run Clomics pipeline."
   echo "Syntax: run_Clomics.sh [-r|p|n|i|d|o|t]"
}

# Get the options
while getopts ":r:p:n:i:d:o:t:" option; do
    case ${option} in
        h) # display Help
            Help
            exit;;
        r) path_to_src=${OPTARG};;
        p) path_to_clomics=${OPTARG};;
        n) n_peaks=${OPTARG};;
        i) input_path=${OPTARG};;
        d) dataset=${OPTARG};;
        o) output_path=${OPTARG};;
        t) thresholds=${OPTARG};;
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done

echo "Running Clomics..."
Rscript \
    ${path_to_src}/main_clomics_script_parallel.R \
    ${path_to_clomics} \
    ${n_peaks} \
    ${dataset} \
    ${input_path}/${dataset}/H3K4me1_${dataset}_RPKMnorm_regressed_qqnorm_forCRD.bed.gz,${input_path}/${dataset}/H3K27ac_${dataset}_RPKMnorm_regressed_qqnorm_forCRD.bed.gz \
    ${output_path} \
    ${thresholds}
