#!/bin/bash
# +
## SPECIFY:

## path to PHM
export phm_dir=/data/software/PHM-0.2;

## path to PHM helper scripts:
path_to_phm_help=/data/pushkare/computational_paper/GITHUB/00.CM_mapping/PHM/PHM_helper;

## output path:
output_path=/data/pushkare/computational_paper/GITHUB/00.CM_mapping/PHM;

## input dataset:
dataset="test_data";

## theshold(s) for posterior probaility
## of causal interaction for DAG construction:
thresholds=( "0.7" "0.75" "0.8" );

## chromosomes to use
chromosomes="chr21,chr22";

## Run PHM on test_data for chr21 and chr22
sh /data/pushkare/computational_paper/GITHUB/00.CM_mapping/PHM/run_PHM.sh \
    -c ${chromosomes} \
    -d ${dataset} \
    -p ${phm_dir} \
    -s ${path_to_phm_help} \
    -t ${thresholds} \
    -o ${output_path}
# -


