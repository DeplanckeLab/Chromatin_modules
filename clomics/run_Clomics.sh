#!/bin/bash
Help()
{
    # Display Help
    echo "Defining Clomics input and running the tool."
    echo "Usage: run_Clomics.sh [-i|p|c|s|d|o]"
    echo "  -i          Path to count matrix/matrices (comma separated)"
    echo "  -p          Path to Clomics tool"
    echo "  -c          Chromosomes to use: 1. string with one chromosome '13'; 2. comma separated chromosomes '1,5,6' (use only chromosomes 1, 5 and 6); 3. dash separated chromosomes '1-22' (all chr in range from 1 to 22)"
    echo "  -s          Path to Clomics helper scripts"
    echo "  -d          Session/dataset name"
    echo "-----------------------------------------------------"
    echo "  -o          Path to output directory"
}

# Get the options
while getopts ":i:p:c:s:d:o:" option; do
    case ${option} in
        h) # display Help
            Help
            exit;;
        i) input_path=${OPTARG};;
        p) path_to_clomics=${OPTARG};;
        c) chromosomes=${OPTARG};;
        s) path_to_scripts=${OPTARG};;
        d) dataset=${OPTARG};;
        o) output_path=${OPTARG};;
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done

echo "Running Clomics..."
Rscript \
    ${path_to_scripts}/main_clomics_script.R \
    ${path_to_clomics} \  ## Path to Clomics tool
    ${input_path} \  ## Normalized matrices in bed.gz format from qtltools (separated by commas)
    ${chromosomes} \  ## Chromosomes to use
    ${output_path}  ## Output folder

echo "Converting CRDs into CM format..."
python3.9 ${path_to_scripts}/CRDs_to_CMs.py \
    ${output_path}/CRD.all.clomics.bed \
    ${output_path}
