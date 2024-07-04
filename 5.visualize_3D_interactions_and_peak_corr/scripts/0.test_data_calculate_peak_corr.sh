dataset="test_data";
path_to_input_directory="/data/pushkare/Chromatin_modules/${dataset}";

k4me1_matrix="H3K4me1:${path_to_input_directory}/H3K4me1_chr22.bed";
k27ac_matrix="H3K27ac:${path_to_input_directory}/H3K27ac_chr22.bed";

path_to_theoretical_pv_output="/data/pushkare/Chromatin_modules/5.visualize_3D_interactions_and_peak_corr/peak_correlations";
scripts_path="/data/pushkare/Chromatin_modules/5.visualize_3D_interactions_and_peak_corr/scripts";

chromosomes="22";  ## "1-22"
n_cores=12;

python3.9 \
    $scripts_path/get_peak_corr.py \
    -d $dataset \
    -f "${k4me1_matrix},${k27ac_matrix}" \
    -c $chromosomes \
    -o $path_to_theoretical_pv_output \
    -n $n_cores

sort -k1,1 -k2,2n \
    ${path_to_theoretical_pv_output}/${dataset}_theoretical_corr_with_p_values.bed \
    > ${path_to_theoretical_pv_output}/${dataset}_theoretical_corr_with_p_values_sorted.bed;
bgzip ${path_to_theoretical_pv_output}/${dataset}_theoretical_corr_with_p_values_sorted.bed;
tabix -p bed ${path_to_theoretical_pv_output}/${dataset}_theoretical_corr_with_p_values_sorted.bed.gz;
rm ${path_to_theoretical_pv_output}/${dataset}_theoretical_corr_with_p_values_sorted.bed;


