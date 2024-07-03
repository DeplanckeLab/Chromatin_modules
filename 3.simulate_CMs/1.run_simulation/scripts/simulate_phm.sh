#!/bin/bash
dataset="test_data";
method="phm";
n_closest_cms=5;
k_sim_per_cm=10;
parallelize=1
n_cores=22;

path_to_src="/data/pushkare/Chromatin_modules/3.simulate_CMs/1.run_simulation/src";
path_to_data="/data/pushkare/Chromatin_modules/3.simulate_CMs";
core_cm_path="/data/pushkare/Chromatin_modules/1.mapped_CMs";

pp_threshold="0.8";

sh ${path_to_src}/simulate.sh \
    -d ${dataset} \
    -m ${method} \
    -i ${path_to_data} \
    -s ${path_to_src} \
    -n ${n_closest_cms} \
    -k ${k_sim_per_cm} \
    -p ${parallelize} \
    -c ${n_cores} \
    -r ${core_cm_path} \
    -t ${pp_threshold}
