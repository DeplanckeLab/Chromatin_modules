#!/bin/bash
# +
## Convert DAGs into VCM-like format

## SPECIFY input dataset
dataset="test_data"

## SPECIFY theshold(s) for posterior probaility
## of causal interaction for DAG construction

thresholds=( "0.7" "0.75" "0.8" );

## SPECIFY path to PHM helper scripts:
path_to_phm_help=/data/pushkare/computational_paper/GITHUB/00.CM_mapping/PHM/PHM_helper;

## SPECIFY OUTPUT PATH
output_path=/data/pushkare/computational_paper/GITHUB/00.CM_mapping/PHM/${dataset}/phm_results;

for threshold in ${thresholds[*]};
do
    python3.9 ${path_to_phm_help}/1.build_DAGs_per_chr.py \
        -d ${dataset} \
        -i ${output_path} \
        -t ${threshold} \
        -c "21-22"
    python3.9 ${path_to_phm_help}/2.get_VCMs_from_DAGs.py \
        -d ${dataset} \
        -i ${output_path} \
        -t ${threshold} \
        -c "21-22"
done
# -


