import os
import numpy as np
import pandas as pd


def check_boundaries(x, y, boundaries_x, boundaries_y):
    if (
        (boundaries_x[0] <= x <= boundaries_x[1])
        and (boundaries_y[0] <= x <= boundaries_y[1])
        and (boundaries_x[0] <= y <= boundaries_x[1])
        and (boundaries_y[0] <= y <= boundaries_y[1])
    ):
        return True
    else:
        return False


# def get_gc_content_per_cm_peak(
#     cm_peaks_dict,
#     gc_content_dict
# ):
#     return {
#         cm_id: [
#             gc_content_dict.get(cm_peak)
#             for cm_peak in cm_peaks
#         ]
#         for cm_id, cm_peaks in cm_peaks_dict.items()
#     }


def subset_tracks_content(
    tracks_df,
    content_df,
    cms_to_subset,
    save_files=True,
    output_path=None,
    output_file_name=None,
):
    tracks_sub = tracks_df.loc[tracks_df.iloc[:, 3].isin(cms_to_subset), :].copy()

    content_sub = content_df.loc[content_df.iloc[:, 0].isin(cms_to_subset), :].copy()
    if save_files:
        tracks_sub.to_csv(
            os.path.join(output_path, output_file_name + ".tracks.bed"),
            sep="\t",
            index=False,
            header=False,
        )
        content_sub.to_csv(
            os.path.join(output_path, output_file_name + ".content.txt"),
            sep="\t",
            index=False,
            header=False,
        )


def get_cm_tracks_content(
    core_cm_path,
    method,
    dataset,
    pp_threshold=None,
    vcm_window=None,
    pv_threshold=None,
    n_peaks=None,
    bg_threshold=None,
):
    if method == "phm":
        method_path = os.path.join(
            core_cm_path,
            method,
            dataset + "_phm_tracks_content",
        )
        tracks_path = os.path.join(
            method_path,
            "_".join([dataset, str(pp_threshold), "merged_phm_all_chr.tracks.bed"]),
        )
        content_path = os.path.join(
            method_path,
            "_".join([dataset, str(pp_threshold), "merged_phm_all_chr.content.txt"]),
        )
    elif method == "vcmtools":
        method_path = os.path.join(core_cm_path, method, "VCMs", vcm_window, dataset)
        tracks_path = os.path.join(
            method_path,
            dataset + "_all_VCMs_corrected_pvalue_" + str(pv_threshold) + ".tracks.bed",
        )
        content_path = os.path.join(
            method_path,
            dataset
            + "_all_VCMs_corrected_pvalue_"
            + str(pv_threshold)
            + ".content.txt",
        )
    elif method == "clomics":
        method_path = os.path.join(core_cm_path, method, n_peaks, bg_threshold, dataset)
        tracks_path = os.path.join(
            method_path, "_".join([dataset, "Clomics_CM.tracks.bed"])
        )
        content_path = os.path.join(
            method_path, "_".join([dataset, "Clomics_CM.content.txt"])
        )
    tracks_df = pd.read_csv(tracks_path, sep="\t", header=None)
    content_df = pd.read_csv(content_path, sep="\t", header=None)

    if not "chr" in content_df.loc[0, 2]:
        for i in range(1, 23):
            content_df.loc[:, 2] = content_df.loc[:, 2].str.replace(
                ":" + str(i) + ":", ":chr" + str(i) + ":"
            )
    return tracks_df, content_df


def get_simulated_cm_tracks_content(
    sim_cm_path, dataset, method, n_closest_cms=5, k_sim_per_cm=10
):
    sim_path = os.path.join(
        sim_cm_path,
        method,
        dataset,
        "_".join(
            [
                str(n_closest_cms),
                "closest_simulated_for",
                str(k_sim_per_cm),
                "sim_per_cm",
            ]
        ),
        "_".join(
            [
                "simulated_" + str(n_closest_cms),
                "closest",
                method,
                dataset,
                "no_repeating_peaks",
            ]
        ),
    )
    sim_tracks_df = pd.read_csv(sim_path + ".tracks.bed", sep="\t", header=None)
    sim_content_df = pd.read_csv(sim_path + ".content.txt", sep="\t", header=None)
    return sim_tracks_df, sim_content_df
