# %%
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy
import sys

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), "src")))
import utils

# %%
dataset = "test_data"
chromosomes = ["chr22"]
data_mapping_dict = {
    dataset: "GM12878",  ## since test_data is from LCLs
}
interaction_type = "MicroC"
resolutions = [500, 1000]

core_path = "/data/pushkare/Chromatin_modules/4.quantify_3D_interactions"

tad_ab_type = "no_TADs_no_AB_compartments"
simulation_type = "stringent_filtering_by_length_and_gc_content"

methods = ["vcmtools", "clomics", "phm"]

output_plots_path = os.path.join(core_path, tad_ab_type, simulation_type, "plots")
if not os.path.exists(output_plots_path):
    os.makedirs(output_plots_path)

# %%
# Define ranges for boxplots
STEP = 25000
ranges = list(
    zip(np.arange(0, 100000, step=STEP), np.arange(STEP, 100000 + STEP, step=STEP))
)
ranges.extend([(100000, 500000)])

# Adjust depending on dataset and resolution:
# format: {"3C_method": {"dataset": {resolution: [y_min, y_max, stars_y_shift]}}}
resolutions_and_y_limits = {
    "HiC": {
        "test_data": {500: [-0.005, 0.05, 0.018], 1000: [-0.5, 7, 0.7]},
    },
    "MicroC": {"test_data": {500: [-0.005, 0.05, 0.01], 1000: [-0.005, 0.08, 0.015]}},
}
# %%
# Get average interactions for reference and simulated CMs (per CM)
for method in methods:
    print(f"Processing {method}...")
    ref_sim_path = os.path.join(
        "/data/pushkare/Chromatin_modules",
        "3.simulate_CMs/",
        "output",
        tad_ab_type,
        "QC",
        simulation_type,
        method,
    )

    for resolution in resolutions:
        if resolution == 500:
            res_str = "500bp"
        elif resolution == 1000:
            res_str = "1kb"
        elif resolution == 5000:
            res_str = "5kb"
        print("Preparing average interactions for", res_str, "resolution")

        input_output_path = os.path.join(
            core_path,
            tad_ab_type,
            simulation_type,
            interaction_type,
            res_str,
            method,
        )

        interactions_by_chr_ref_cm = utils.load_interactions_per_chr(
            input_output_path,
            interaction_type,
            "reference",
            dataset,
            method,
            data_mapping_dict,
            resolution,
        )

        interactions_by_chr_sim_cm = utils.load_interactions_per_chr(
            input_output_path,
            interaction_type,
            "simulated",
            dataset,
            method,
            data_mapping_dict,
            resolution,
        )

        sim_cm_tracks = utils.read_tracks_df(
            ref_sim_path,
            dataset,
            method,
            "sim",
            file_suffix="length_GC_content_filtering",
        )

        if "chr" not in str(sim_cm_tracks.loc[0, 0]):
            sim_cm_tracks.loc[:, 0] = "chr" + sim_cm_tracks.loc[:, 0].astype(str)

        chr_ref_sim_interaction_dict = {}
        for chromosome in chromosomes:
            chr_df = sim_cm_tracks.loc[sim_cm_tracks.loc[:, 0] == chromosome, :]
            chr_dict = {
                ref_cm: list(ref_sim_df.loc[:, 3])
                for ref_cm, ref_sim_df in chr_df.loc[:, [13, 3]].groupby(13)
            }
            ref_sim_interaction_freq = {}
            for ref_param_dict in interactions_by_chr_ref_cm.get(chromosome):
                sim_values = [
                    sim_param_dict["avg_hic_frequency"]
                    for sim_cm in chr_dict.get(ref_param_dict["cm_id"])
                    for sim_param_dict in interactions_by_chr_sim_cm.get(chromosome)
                    if sim_cm in sim_param_dict["cm_id"]
                ]
                ref_sim_interaction_freq[ref_param_dict["cm_id"]] = {
                    "ref_cm_length": ref_param_dict["cm_length"],
                    "ref_hic_frequency": ref_param_dict["avg_hic_frequency"],
                    "avg_sim_hic_frequency": np.mean(sim_values),
                }
            chr_ref_sim_interaction_dict[chromosome] = ref_sim_interaction_freq
        utils.save_average_interactions(
            chr_ref_sim_interaction_dict,
            input_output_path,
            dataset,
            method,
            data_mapping_dict,
            interaction_type,
            resolution,
            hichip_type=None,
        )
print("DONE!")
# %%
# Get average interactions by ranges and plot boxplots
for method in methods:
    for resolution in resolutions:
        y_min, y_max, stars_y_shift = resolutions_and_y_limits[interaction_type][
            dataset
        ][resolution]
        if resolution == 500:
            res_str = "500bp"
        elif resolution == 1000:
            res_str = "1kb"
        elif resolution == 5000:
            res_str = "5kb"
        input_output_path = os.path.join(
            core_path, tad_ab_type, simulation_type, interaction_type, res_str, method
        )

        chr_ref_sim_interaction_dict = utils.load_average_interactions(
            input_output_path,
            dataset,
            method,
            data_mapping_dict,
            interaction_type,
            resolution,
            hichip_type=None,
        )
        ds_ranges = utils.prepare_data_by_ranges(chr_ref_sim_interaction_dict, ranges)
        utils.boxplot_by_range(
            ds_ranges,
            dataset,
            method,
            interaction_type,
            res_str,
            vicinity=False,
            save_figure=True,
            path=output_plots_path,
            file_name="_".join(
                [interaction_type, dataset, method, "by_ranges", res_str + ".pdf"]
            ),
            y_min=y_min,
            y_max=y_max,
            stars_y_shift=stars_y_shift,
        )
# %%
