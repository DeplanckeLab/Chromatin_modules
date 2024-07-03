Help()
{
   ## Display Help
   echo "Run CM simulation pipeline."
   echo "Syntax: simulate.sh [-d|m|i|s|n|k|p|c|r|t|w|v|a|b|h]"
}

n_peaks="None";
bg_threshold="None";
pp_threshold="None";
vcm_window="None";
pv_threshold="None";

## Get the options
while getopts ":d:m:i:s:n:k:p:c:r:t:w:v:a:b:" option; do
    case $option in
        h) ## display Help
            Help
            exit;;
        d) dataset=${OPTARG};;
        m) method=${OPTARG};;
        i) path_to_data=${OPTARG};;
        s) path_to_src=${OPTARG};;
        n) n_closest_cms=${OPTARG};;
        k) k_sim_per_cm=${OPTARG};;
        p) parallelize=${OPTARG};;
        c) n_cores=${OPTARG};;
        r) core_cm_path=${OPTARG};;
        t) pp_threshold=${OPTARG};;
        w) vcm_window=${OPTARG};;
        v) pv_threshold=${OPTARG};;
        a) n_peaks=${OPTARG};;
        b) bg_threshold=${OPTARG};;
        \?) ## Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done

path_to_input_directory=${path_to_data}/input/${method};
path_to_output_directory=${path_to_data}/output/no_TADs_no_AB_compartments;

mkdir -p ${path_to_output_directory}/${method}/${dataset}/${n_closest_cms}_closest_simulated_for_${k_sim_per_cm}_sim_per_cm;
python3 ${path_to_src}/0.simulate_CMs.py \
    ${dataset} \
    ${path_to_input_directory}/${dataset} \
    ${path_to_output_directory}/${method}/${dataset}/${n_closest_cms}_closest_simulated_for_${k_sim_per_cm}_sim_per_cm \
    ${method} \
    ${n_closest_cms} \
    ${k_sim_per_cm} \
    ${parallelize} \
    ${n_cores}

python3 ${path_to_src}/1.filter_simulated_CMs.py \
    -d ${dataset} \
    -m ${method} \
    --input_path ${path_to_output_directory} \
    --core_cm_path ${core_cm_path} \
    --n_closest_cms ${n_closest_cms} \
    --k_sim_per_cm ${k_sim_per_cm} \
    --pp_threshold ${pp_threshold} \
    --vcm_window ${vcm_window} \
    --pv_threshold ${pv_threshold} \
    --n_peaks ${n_peaks} \
    --bg_threshold ${bg_threshold}