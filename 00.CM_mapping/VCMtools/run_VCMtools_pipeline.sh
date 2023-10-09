dataset="test_data";
path_to_input_directory=/data/pushkare/computational_paper/GITHUB/test_data;
path_to_theoretical_pv_output=/data/pushkare/computational_paper/GITHUB/00.CM_mapping/VCMtools/output/corr_files/theoretical_corr_output;
path_to_empirical_pv_output=/data/pushkare/computational_paper/GITHUB/00.CM_mapping/VCMtools/output/corr_files/empirical_corr_output;
scripts_path=/data/pushkare/computational_paper/GITHUB/00.CM_mapping/VCMtools/main;
vcm_output_path=/data/pushkare/computational_paper/GITHUB/00.CM_mapping/VCMtools/output/CMs;
chromosomes="21,22";
n_cores=2;

python3.9 \
    $scripts_path/01.VCMtools_theoretical_pv.py \
    -d $dataset \
    -i $path_to_input_directory \
    -c $chromosomes \
    -o $path_to_theoretical_pv_output \
    -n $n_cores

python3.9 \
    $scripts_path/02.VCMtools_empirical_pv.py \
    -d $dataset \
    -i $path_to_theoretical_pv_output \
    -o $path_to_empirical_pv_output \

python3.9 \
    $scripts_path/03.get_VCMs_from_edgelist.py \
    -d $dataset \
    -corr_p $path_to_empirical_pv_output/${dataset}_empirical_corr_p_values.txt \
    -pv 0.001 \
    -t 1 \
    -f 1 \
    -m 1 \
    -o $vcm_output_path

