#!/bin/bash
Help()
{
   # Display Help
   echo "Run VCMtools python pipeline."
   echo "Syntax: pyVCMtools.sh [-s|d|m|p|f|c||n|o]"
}

# Get the options
while getopts ":h:s:d:m:p:f:c:n:o:" option; do
    case $option in
        h) # display Help
            Help
            exit;;
        s) scripts_path=${OPTARG};;
        d) dataset=${OPTARG};;
        m) max_peak_dist=${OPTARG};;
        p) pv_thresholds=${OPTARG};;
        f) input_files=${OPTARG};;
        c) chromosomes=${OPTARG};;
        n) n_cores=${OPTARG};;
        o) output_path=${OPTARG};;
        
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done

theoretical_corr_path=${output_path}/theoretical_corr_pv/${max_peak_dist}Mb;
empirical_corr_path=${output_path}/empirical_corr_pv/${max_peak_dist}Mb;
vcm_path=${output_path}/VCMs/${max_peak_dist}Mb/${dataset};
mkdir -p ${theoretical_corr_path};
mkdir -p ${empirical_corr_path};
mkdir -p ${vcm_path};

python3 \
    ${scripts_path}/01.VCMtools_theoretical_pv.py \
    --dataset $dataset \
    --max_peak_dist $max_peak_dist \
    --input_files $input_files \
    --chromosomes $chromosomes \
    --output_path ${theoretical_corr_path} \
    --n_cores $n_cores

python3 \
    ${scripts_path}/02.VCMtools_empirical_pv.py \
    --dataset $dataset \
    --path_to_input_directory ${theoretical_corr_path} \
    --output_path ${empirical_corr_path}

python3 \
    ${scripts_path}/03.get_VCMs_from_edgelist.py \
    --dataset $dataset \
    --correlation_path ${empirical_corr_path}/${dataset}_all_marks_empirical_corr_p_values.txt \
    --pvalue ${pv_thresholds} \
    --output_path ${vcm_path}
