#!/bin/bash
Help()
{
   # Display Help
   echo "Run VCMtools python pipeline."
   echo "Syntax: VCMtools.sh [-s|d|f|c|p|n|o]"
}

# +
pvalue=0.001;
n_cores=1;

# Get the options
while getopts ":h:s:d:f:c:p:n:o:" option; do
    case $option in
        h) # display Help
            Help
            exit;;
        s) scripts_path=${OPTARG};;
        d) dataset=${OPTARG};;
        f) path_to_input=${OPTARG};;
        c) chromosomes=${OPTARG};;
        p) pvalue=${OPTARG};;
        n) n_cores=${OPTARG};;
        o) output_path=${OPTARG};;
        
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done

# +
mkdir -p ${output_path}/corr_files/theoretical_pv_output
mkdir -p ${output_path}/corr_files/empirical_pv_output
mkdir -p ${output_path}/CMs

python3.9 \
    ${scripts_path}/01.get_corr_and_theoretical_pv.py \
    -d $dataset \
    -i $path_to_input \
    -c $chromosomes \
    -o ${output_path}/corr_files/theoretical_pv_output \
    -n $n_cores

python3.9 \
    ${scripts_path}/02.get_empirical_pv.py \
    -d $dataset \
    -i ${output_path}/corr_files/theoretical_pv_output \
    -o ${output_path}/corr_files/empirical_pv_output \

python3.9 \
    ${scripts_path}/03.get_CMs_from_edgelist.py \
    -d $dataset \
    -corr_p ${output_path}/corr_files/empirical_pv_output/${dataset}_empirical_corr_p_values.txt \
    -pv ${pvalue} \
    -t 1 \
    -f 1 \
    -m 1 \
    -o $output_path/CMs
