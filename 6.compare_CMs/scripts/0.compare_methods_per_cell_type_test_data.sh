#!/bin/bash
dataset="test_data"
score_type="HM";

## Clomics parameters
n_peaks="n_peaks_200";
bg_threshold="bg_threshold_3";
## VCMtools parameters
vcm_window="0.5Mb";
pv_threshold="0.001";
## PHM parameter
pp_threshold="0.8";

## Use a range of similarity thresholds
similarity_threshold="0-1;0.1";

## Input paths
cm_path="/data/pushkare/Chromatin_modules/1.mapped_CMs";
cm_peak_path="/data/pushkare/Chromatin_modules/2.peaks_in_CMs/peak_files";
core_path=/data/pushkare/Chromatin_modules/6.compare_CMs;
src_path=${core_path}/src;

compare="methods";
if [[ ${compare} == "methods" ]]
then
    output_path=${core_path}/compare_methods_per_cell_type;
else
    output_path=${core_path}/compare_cell_types_per_method;
fi

ds_output_path=${output_path}/${dataset};
score_output_path=${ds_output_path}/${score_type};
mkdir -p ${score_output_path};

# Iterate over all possible method pairs
for i in "vcmtools clomics" "vcmtools phm" "clomics vcmtools" "clomics phm" "phm vcmtools" "phm clomics" "vcmtools vcmtools" "clomics clomics" "phm phm";
do
    method_pair=( $i );
    method_1=${method_pair[0]}
    method_2=${method_pair[1]}

    if [[ ${method_1} == "clomics" ]]
    then
        suffix_1="_all";
    else
        suffix_1="_NOT_TOTEM";
    fi

    if [[ ${method_2} == "clomics" ]]
    then
        suffix_2="_all";
    else
        suffix_2="_NOT_TOTEM";
    fi
    # Overlap peaks in chromatin modules from different cell types with each other
    if ! [ -f "${ds_output_path}/${method_1}_${method_2}_jaccard_peak_overlap.bed" ]
        then
            echo "Preparing the peak-based chromatin module ovelap files..."
            file_1=${cm_peak_path}/${dataset}_${method_1}${suffix_1}_peaks.bed
            file_2=${cm_peak_path}/${dataset}_${method_2}${suffix_2}_peaks.bed
            bedtools intersect -wo \
                -a ${file_1} \
                -b ${file_2} \
                > ${ds_output_path}/${method_1}_${method_2}_jaccard_peak_overlap.bed
    fi
done

python3 ${src_path}/main_cm_overlap_per_cell_type.py \
    -m ${dataset} \
    -d "vcmtools,clomics,phm" \
    -cm ${cm_path} \
    -pcm ${cm_peak_path} \
    -pp ${pp_threshold} \
    -vw ${vcm_window} \
    -pv ${pv_threshold} \
    -np ${n_peaks} \
    -bg ${bg_threshold} \
    -pop ${ds_output_path} \
    --use_cmQTLs 0 \
    --min_cm_size 0 \
    -s ${score_type} \
    --similarity_threshold ${similarity_threshold} \
    -o ${score_output_path} \
    --overwrite 1

