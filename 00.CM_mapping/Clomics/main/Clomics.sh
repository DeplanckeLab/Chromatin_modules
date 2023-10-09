#!/bin/bash
Help()
{
   # Display Help
   echo "Run Clomics pipeline."
   echo "Syntax: run_Clomics.sh [-i|s|d|m|o]"
}

marks="all_marks";

# Get the options
while getopts ":i:s:d:m:o:" option; do
    case ${option} in
        h) # display Help
            Help
            exit;;
        i) input_path=${OPTARG};;  # Path to count matrix/matrices (comma separated)
        s) path_to_scripts=${OPTARG};;
        d) dataset=${OPTARG};;
        m) marks=${OPTARG};;
        o) output_path=${OPTARG};;
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done

clomics_out_path=${output_path}/Clomics/${marks}/${dataset};
mkdir -p ${output_path}

sh ${path_to_scripts}/run_Clomics.sh \
    -i ${input_path} \
    -d ${dataset} \
    -o ${clomics_out_path}

echo "Running Clomics..."
Rscript \
    ${path_to_scripts}/main_clomics_script.R \
    ${input_path} \
    ${clomics_out_path}

echo "Converting CRDs into VCM-like format..."
python3.9 ${path_to_scripts}/CRDs_to_CMs.py \
    ${clomics_out_path}/CRD.all.clomics.bed \
    ${clomics_out_path}
