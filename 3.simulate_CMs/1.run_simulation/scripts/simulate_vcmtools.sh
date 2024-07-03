#!/bin/bash
dataset="test_data";
method="vcmtools";
n_closest_cms=5;
k_sim_per_cm=10;
parallelize=1
n_cores=22;

core_cm_path="/data/pushkare/Chromatin_modules/1.mapped_CMs";
path_to_data="/data/pushkare/Chromatin_modules/3.simulate_CMs";
path_to_src="/data/pushkare/Chromatin_modules/3.simulate_CMs/1.run_simulation/src";

vcm_window="0.5Mb";
pv_threshold="0.001";

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
    -w ${vcm_window} \
    -v ${pv_threshold}
