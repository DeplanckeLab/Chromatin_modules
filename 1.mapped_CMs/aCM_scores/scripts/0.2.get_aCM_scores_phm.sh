#!/bin/bash
dataset="test_data";

## Specify method and parameters
method="phm";
pp_threshold="0.8";
method_specific_path="${dataset}_phm_tracks_content";
method_specific_ext="0.8_merged_phm_all_chr";

core_path="/data/pushkare/Chromatin_modules/1.mapped_CMs";
path_to_scripts=${core_path}/aCM_scores/scripts/src;

path_to_data="/data/pushkare/Chromatin_modules/test_data";
path_to_output=${core_path}/aCM_scores/${method}/${dataset};

mkdir -p ${path_to_output};
mkdir -p ${path_to_output}/plot
mkdir -p ${path_to_output}/aCM_matrix;

Rscript ${path_to_scripts}/calculate_aCM.R \
    ${core_path}/${method}/${method_specific_path}/${dataset}_${method_specific_ext}.content.txt \
    ${core_path}/${method}/${method_specific_path}/${dataset}_${method_specific_ext}.tracks.bed \
    "H3K4me1:${path_to_data}/gzipped_bed/H3K4me1_chr22.bed.gz,H3K27ac:${path_to_data}/gzipped_bed/H3K27ac_chr22.bed.gz" \
    ${path_to_data}/LCL_genotypes_chr22.vcf.gz \
    ${path_to_output}/${method}_${dataset}_merged_count_matrix.bed \
    ${path_to_output}/plot \
    ${path_to_output}/aCM_matrix \
    ${path_to_output}/${method}_${dataset}_cmQTLs.txt \
    ${dataset}

## If Error in smooth.spline(lambda, pi0, df = smooth.df): 
## missing or infinite values in inputs are not allowed in ->
## not enough p-values! can be fixed by changing line 90 in qtltools_runFDR_cis.R
## to Q <- qvalue(D[,opt_col], pi0 = 1)
Rscript ${path_to_scripts}/qtltools_runFDR_cis.R \
    ${path_to_output}/${method}_${dataset}_cmQTLs.txt \
    0.05 \
    ${path_to_output}/${method}_${dataset}_cmQTLs
