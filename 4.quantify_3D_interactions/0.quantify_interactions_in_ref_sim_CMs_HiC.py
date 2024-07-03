# %%
import cooler
import numpy as np
import os
import pandas as pd
import sys
import tqdm.auto as tqdm

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), "src")))
from HiCM import HiCM
import utils

# %% [markdown]
# ### Data preparation

# %%
chromosomes = ["chr22"]

dataset = "test_data"
methods = ["vcmtools", "clomics", "phm"]
save_files = True

hic_data = "GM12878"
interaction_type = "MicroC"  # "HiC"
resolutions = [500, 1000]

chr_sizes_path = (
    "/data/pushkare/Chromatin_modules/genome_annotations/chr_sizes_hg19.txt"
)
core_path = "/data/pushkare/Chromatin_modules/4.quantify_3D_interactions"

tad_ab_type = "no_TADs_no_AB_compartments"
simulation_type = "stringent_filtering_by_length_and_gc_content"
ref_sim_file_suffix = "length_GC_content_filtering"

path_to_3d_data = os.path.join(
    "/data/pushkare/computational_paper",
    "08.3D_interactions_of_CREs_in_chromatin_modules",
    "3D_chromatin_data",
    interaction_type,
)

if interaction_type == "HiC":
    cool_dict = {
        "500bp": path_to_3d_data + "_500bp/LCL_mega_42B_500bp_30_cool.cool",
        "1kb": os.path.join(
            path_to_3d_data,
            "Rao2014-" + hic_data + "-MboI-allreps-filtered.1kb.cool",
        ),
    }
elif interaction_type == "MicroC":
    cool_dict = {
        "500bp": os.path.join(path_to_3d_data, "microc_800m_500bp.cool"),
        "1kb": os.path.join(path_to_3d_data, "microc_800m_1kb.cool"),
    }

# %%
for method in methods:
    print("Processing", method, "...")
    ref_sim_path = os.path.join(
        "/data/pushkare/Chromatin_modules",
        "3.simulate_CMs",
        "output",
        tad_ab_type,
        "QC",
        simulation_type,
        method,
    )
    ref_files_dict = {
        "tracks": os.path.join(
            ref_sim_path,
            "_".join(
                [
                    dataset,
                    "ref",
                    method,
                    ref_sim_file_suffix + ".tracks.bed",
                ]
            ),
        ),
        "content": os.path.join(
            ref_sim_path,
            "_".join(
                [
                    dataset,
                    "ref",
                    method,
                    ref_sim_file_suffix + ".content.txt",
                ]
            ),
        ),
    }
    sim_files_dict = {
        "tracks": os.path.join(
            ref_sim_path,
            "_".join(
                [
                    dataset,
                    "sim",
                    method,
                    ref_sim_file_suffix + ".tracks.bed",
                ]
            ),
        ),
        "content": os.path.join(
            ref_sim_path,
            "_".join(
                [
                    dataset,
                    "sim",
                    method,
                    ref_sim_file_suffix + ".content.txt",
                ]
            ),
        ),
    }
    for resolution in resolutions:
        if resolution == 500:
            res_str = "500bp"
        elif resolution == 1000:
            res_str = "1kb"
        elif resolution == 5000:
            res_str = "5kb"

        print("Preparing interactions for", res_str, "resolution")
        if not os.path.exists(
            os.path.join(
                core_path,
                tad_ab_type,
                simulation_type,
                interaction_type,
                res_str,
                method,
            )
        ):
            os.makedirs(
                os.path.join(
                    core_path,
                    tad_ab_type,
                    simulation_type,
                    interaction_type,
                    res_str,
                    method,
                )
            )

        cool_file = cool_dict.get(res_str)
        cool_mtx = cooler.Cooler(cool_file)
        hiCM = HiCM(
            cool_mtx=cool_mtx, resolution=resolution, chr_sizes_path=chr_sizes_path
        )

        ref_cm_data = utils.get_content_dict(
            ref_files_dict.get("content"),
        )

        sim_cm_data = utils.get_content_dict(
            sim_files_dict.get("content"),
        )
        peaks_per_ref_cm = ref_cm_data["cm_content"]
        peaks_per_sim_cm = sim_cm_data["cm_content"]

        interactions_by_chr_ref_cm = utils.get_hic_interactions(
            hiCM, peaks_per_ref_cm, cool_mtx.chromnames[:-3]
        )
        interactions_by_chr_sim_cm = utils.get_hic_interactions(
            hiCM, peaks_per_sim_cm, cool_mtx.chromnames[:-3]
        )
        if save_files:
            np.save(
                os.path.join(
                    core_path,
                    tad_ab_type,
                    simulation_type,
                    interaction_type,
                    res_str,
                    method,
                    "_".join(
                        [
                            dataset,
                            "reference",
                            method,
                            hic_data,
                            interaction_type,
                            "interactions",
                            str(resolution) + "bp.npy",
                        ]
                    ),
                ),
                interactions_by_chr_ref_cm,
                allow_pickle=True,
            )
            np.save(
                os.path.join(
                    core_path,
                    tad_ab_type,
                    simulation_type,
                    interaction_type,
                    res_str,
                    method,
                    "_".join(
                        [
                            dataset,
                            "simulated",
                            method,
                            hic_data,
                            interaction_type,
                            "interactions",
                            str(resolution) + "bp.npy",
                        ]
                    ),
                ),
                interactions_by_chr_sim_cm,
                allow_pickle=True,
            )
print("DONE!")
# %%
