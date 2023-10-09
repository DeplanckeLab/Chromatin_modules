#!/bin/bash
Help()
{
   # Display Help
   echo "Run VCMtools python pipeline."
   echo "Syntax: pyVCMtools.sh [-s|d|f|c|t|e|n|o]"
}

# Get the options
while getopts ":h:s:d:f:c:t:e:n:o:" option; do
    case $option in
        h) # display Help
            Help
            exit;;
        s) scripts_path=${OPTARG};;
        d) dataset=${OPTARG};;
        f) input_files=${OPTARG};;
        c) chromosomes=${OPTARG};;
        t) path_to_theoretical_pv_output=${OPTARG};;
        e) path_to_empirical_pv_output=${OPTARG};;
        n) n_cores=${OPTARG};;
        o) vcm_output_path=${OPTARG};;
        
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
   esac
done

python3.9 \
    ${scripts_path}/01.VCMtools_theoretical_pv.py \
    -d $dataset \
    -f $input_files \
    -c $chromosomes \
    -o $path_to_theoretical_pv_output \
    -n $n_cores

python3.9 \
    ${scripts_path}/02.VCMtools_empirical_pv.py \
    -d $dataset \
    -i $path_to_theoretical_pv_output \
    -o $path_to_empirical_pv_output \

python3.9 \
    ${scripts_path}/03.get_VCMs_from_edgelist.py \
    -d $dataset \
    -corr_p ${path_to_empirical_pv_output}/${dataset}_all_marks_empirical_corr_p_values.txt \
    -pv 0.001 \
    -t 1 \
    -f 1 \
    -m 1 \
    -o $vcm_output_path
