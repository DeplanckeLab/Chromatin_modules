import pandas as pd
import os
import sys
import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib

import tqdm.auto as tqdm


def get_content_dict(
    content_file,
    #     tracks_file
):
    # Upload .content file to get CM peaks
    cm_content = pd.read_csv(
        content_file, sep="\t", header=None, names=["cm_id", "cm_size", "peaks"]
    )
    if not "chr" in cm_content.iloc[0, 2]:
        for i in range(1, 23):
            cm_content.loc[:, "peaks"] = cm_content.loc[:, "peaks"].str.replace(
                ":" + str(i) + ":", ":chr" + str(i) + ":"
            )
    #     # Upload .content file to get CM peaks
    #     cm_tracks = pd.read_csv(tracks_file, sep='\t', header=None)
    #     if not 'chr' in str(cm_tracks.iloc[0, 0]):
    #         cm_tracks.loc[:, 0] = 'chr' + cm_tracks.loc[:, 0].astype(str)

    cm_content_dict = dict(
        zip(cm_content.loc[:, "cm_id"], cm_content.loc[:, "peaks"].str.split(","))
    )
    return {
        "cm_content": {
            cm_id: [":".join(cm_peak.split(":")[-3:]) for cm_peak in cm_peaks]
            for cm_id, cm_peaks in cm_content_dict.items()
        }
    }


def get_hic_interactions(hicm_object, peaks_dict, chromosomes):

    interactions_by_chr_cm = {}
    for chromosome in chromosomes:
        if "chr" not in chromosome:
            chromosome = "chr" + chromosome
        interactions_by_chr_cm[chromosome] = []
    for cm_id, cm_peaks in tqdm.tqdm(peaks_dict.items()):
        chromosome_id = cm_peaks[0].split(":")[-3]
        if "chr" in chromosome_id:
            chrom_id = chromosome_id.split("chr")[1]
        else:
            chrom_id = chromosome_id

        cm_peak_coord = sorted(
            [[int(peak.split(":")[-2]), int(peak.split(":")[-1])] for peak in cm_peaks]
        )
        cm_start = min(cm_peak_coord, key=lambda x: x[0])[0]
        cm_end = max(cm_peak_coord, key=lambda x: x[1])[1]
        avg_frequencies = hicm_object.get_peak_interactions(
            chrom_id,
            cm_peak_coord,
            drop_zeros=True,
        )

        interactions_by_chr_cm[chromosome_id].append(
            {
                "cm_id": cm_id,
                "avg_hic_frequency": avg_frequencies,
                "cm_length": cm_end - cm_start,
            }
        )
    return interactions_by_chr_cm


def merge_intersecting_intervals(list_of_pairs):
    """For the given list of pairs the function finds and merges intersecting intervals.
    The output is the list of intervals of length <= length(input list)"""
    list_of_pairs = sorted(list_of_pairs)
    merged_intervals = []
    for pair in list_of_pairs:
        if not merged_intervals:
            merged_intervals.append((pair[0], pair[1]))
        else:
            lower_pair = merged_intervals[-1]
            # test for intersection between lower and higher bounds
            if pair[0] <= lower_pair[1]:
                upper_bound = max(lower_pair[1], pair[1])
                # replace with the merged interval
                merged_intervals[-1] = (lower_pair[0], upper_bound)
            else:
                merged_intervals.append((pair[0], pair[1]))
    return merged_intervals


def get_values_by_dataset(dataset_dict, a, b, vicinity=False, c=None, d=None):
    values_by_dataset = []
    significant_differences = []
    stars = []
    for pos, (hichip_mark, ref_sim_dict) in enumerate(dataset_dict.items()):
        ref_lst = np.log2(np.array(ref_sim_dict["reference"]) + 1)
        sim_lst = np.log2(np.array(ref_sim_dict["simulated"]) + 1)
        if vicinity:
            ref_vicinity_lst = np.log2(np.array(ref_sim_dict["reference vicinity"]) + 1)
            sim_vicinity_lst = np.log2(np.array(ref_sim_dict["simulated vicinity"]) + 1)
            values_by_dataset.append(ref_vicinity_lst[~np.isnan(ref_vicinity_lst)])
        mask = ~np.isnan(ref_lst) & ~np.isnan(sim_lst)
        ref_lst = ref_lst[mask]
        sim_lst = sim_lst[mask]
        values_by_dataset.append(ref_lst)
        values_by_dataset.append(sim_lst)

        if vicinity:
            values_by_dataset.append(sim_vicinity_lst[~np.isnan(sim_vicinity_lst)])
        if len(ref_lst) < 2:
            ref_sim_pv = 1
        else:
            if len(ref_lst) != len(sim_lst):
                print(
                    "Lengths of reference and simulated interaction lists are different! Exiting..."
                )
                sys.exit()
            if all(ref_lst - sim_lst) == 0:
                ref_sim_pv = 1
            else:
                _, ref_sim_pv = scipy.stats.wilcoxon(
                    ref_lst, sim_lst, correction=True, alternative="greater"
                )
        if vicinity:
            _, ref_v_sim_pv = scipy.stats.mannwhitneyu(
                ref_vicinity_lst[~np.isnan(ref_vicinity_lst)],
                sim_lst[~np.isnan(sim_lst)],
            )
            _, ref_sim_v_pv = scipy.stats.mannwhitneyu(
                ref_lst[~np.isnan(ref_lst)],
                sim_vicinity_lst[~np.isnan(sim_vicinity_lst)],
            )
            _, ref_v_sim_v_pv = scipy.stats.mannwhitneyu(
                ref_vicinity_lst[~np.isnan(ref_vicinity_lst)],
                sim_vicinity_lst[~np.isnan(sim_vicinity_lst)],
            )
            _, ref_ref_v_pv = scipy.stats.mannwhitneyu(
                ref_lst[~np.isnan(ref_lst)],
                ref_vicinity_lst[~np.isnan(ref_vicinity_lst)],
            )
            _, sim_sim_v_pv = scipy.stats.mannwhitneyu(
                sim_lst[~np.isnan(sim_lst)],
                sim_vicinity_lst[~np.isnan(sim_vicinity_lst)],
            )
            if ref_sim_pv < 0.05:
                significant_differences.append([a[pos], c[pos]])
            if ref_v_sim_pv < 0.05:
                significant_differences.append([b[pos], c[pos]])
            if ref_sim_v_pv < 0.05:
                significant_differences.append([a[pos], d[pos]])
            if ref_v_sim_v_pv < 0.05:
                significant_differences.append([b[pos], d[pos]])
            if ref_ref_v_pv < 0.05:
                significant_differences.append([a[pos], b[pos]])
            if sim_sim_v_pv < 0.05:
                significant_differences.append([c[pos], d[pos]])
        else:
            if ref_sim_pv <= 0.05:
                significant_differences.append([a[pos], b[pos]])
                if ref_sim_pv <= 0.001:
                    stars.append("***")
                elif (ref_sim_pv > 0.001) and (ref_sim_pv <= 0.01):
                    stars.append("**")
                elif (ref_sim_pv > 0.01) and (ref_sim_pv <= 0.05):
                    stars.append("*")
    return values_by_dataset, significant_differences, stars


def boxplot_by_range(
    chr_dict,
    dataset,
    method,
    interaction_type,
    resolution,
    vicinity=False,
    save_figure=False,
    path=None,
    file_name=None,
    y_min=None,
    y_max=None,
    stars_y_shift=None,
    hichip_type=None,
):
    labels = [
        str(
            str(round(length_range[0] / 1000))
            + "-"
            + str(round(length_range[1] / 1000))
            + "kb"
        )
        for length_range in chr_dict.keys()
    ]
    a = np.arange(1, len(labels) * 2 + 1, 2)
    b = a + 0.4
    if vicinity:
        c = b + 0.25
        d = c + 0.25
    else:
        c, d = None, None

    positions = []
    if vicinity:
        x_ticks = b + 0.125
    else:
        x_ticks = a + 0.125

    for i, el_a in enumerate(a):
        positions.append(el_a)
        positions.append(b[i])
        if vicinity:
            positions.append(c[i])
            positions.append(d[i])

    values_by_dataset, significant_differences, stars = get_values_by_dataset(
        chr_dict, a, b, vicinity=False, c=None, d=None
    )
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    bplot = ax.boxplot(
        values_by_dataset,
        positions=positions,
        widths=0.3,
        #         notch=True,
        showfliers=False,
        vert=True,
        patch_artist=True,
        flierprops=dict(marker=".", markersize=3),
        medianprops=dict(color="black"),
        #         showmeans=True,
        #         meanprops=dict(
        #             markerfacecolor='crimson',
        #             markeredgecolor='crimson',
        #             marker='s',
        #             markersize=3
        #         )
    )

    trans = matplotlib.transforms.blended_transform_factory(ax.transData, ax.transData)
    for s, (start, end) in enumerate(significant_differences):
        max_1 = np.quantile(values_by_dataset[round(start)], 0.75) + 1.5 * (
            np.quantile(values_by_dataset[round(start)], 0.75)
            - np.quantile(values_by_dataset[round(start)], 0.25)
        )
        sorted_lst_1 = np.array(sorted(values_by_dataset[round(start)]))[
            sorted(values_by_dataset[round(start)]) <= max_1
        ]
        whisker_1 = sorted_lst_1[len(sorted_lst_1) - 1]

        max_2 = np.quantile(values_by_dataset[round(end)], 0.75) + 1.5 * (
            np.quantile(values_by_dataset[round(end)], 0.75)
            - np.quantile(values_by_dataset[round(end)], 0.25)
        )

        sorted_lst_2 = np.array(sorted(values_by_dataset[round(end)]))[
            sorted(values_by_dataset[round(end)]) <= max_2
        ]
        whisker_2 = sorted_lst_2[len(sorted_lst_2) - 1]

        q75 = max(whisker_1, whisker_2) + stars_y_shift
        ax.annotate(
            stars[s],
            xy=(start + 0.23, q75),
            xytext=(start + 0.23, q75),
            transform=trans,
            ha="center",
            va="bottom",
            arrowprops=dict(
                arrowstyle="-[, widthB=" + str((end - start) * 1.8) + ", lengthB=0.5",
                lw=0.7,
            ),
        )
    data_name = dataset
    if interaction_type == "HiChIP":
        ax.set_title(
            " ".join(
                [
                    interaction_type,
                    hichip_type,
                    "interaction frequencies for reference and simulated CMs by range,\n",
                    data_name + ", " + method + ",",
                    resolution + "\n",
                ]
            ),
            size=12,
        )
    else:
        ax.set_title(
            " ".join(
                [
                    interaction_type,
                    "interaction frequencies for reference and simulated CMs by range,\n",
                    data_name + ", " + method + ",",
                    resolution + "\n",
                ]
            ),
            size=12,
        )
    if vicinity:
        color_lst = ["lightgray", "lightblue", "darkorange", "gold"]
    else:
        color_lst = ["whitesmoke", "lightgray"]
    colors = color_lst * len(chr_dict.keys())
    for i, (patch, color) in enumerate(zip(bplot["boxes"], colors)):
        #         patch.set_facecolor(color)
        if i % 2 == 0:
            patch.set(
                facecolor=color, edgecolor="black", linewidth=1, alpha=1, hatch="///"
            )
        else:
            patch.set(
                facecolor=color,
                edgecolor="black",
                linewidth=1,
                alpha=1,
            )

    if vicinity:
        ax.legend(
            [
                bplot["boxes"][0],
                bplot["boxes"][1],
                bplot["boxes"][2],
                bplot["boxes"][3],
            ],
            ["Reference", "Reference vicinity", "Simulated", "Simulated vicinity"],
            bbox_to_anchor=(1.5, 0.5),
        )
    else:
        ax.legend(
            [bplot["boxes"][0], bplot["boxes"][1]],
            ["Reference", "Simulated"],
            bbox_to_anchor=(1.5, 0.5),
        )
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(labels, size=12, rotation=45)
    ax.set_xlabel("CM length intervals", size=12)
    ax.set_ylabel("$\\log_2$(mean " + interaction_type + " interaction)", size=12)
    ax.set_ylim((y_min, y_max))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    for i, (position, values_in_region) in enumerate(zip(positions, values_by_dataset)):
        if i % 2 == 0:
            n_points = len(values_in_region)
            plt.text(
                position + 0.2,
                y_max,
                f"n={n_points}",
                ha="center",
                va="bottom",
                size=10,
            )

    if save_figure:
        if not os.path.exists(path):
            os.makedirs(path)
        plt.savefig(
            os.path.join(path, file_name),
            dpi=300,
            transparent=True,
            bbox_inches="tight",
        )


def save_average_interactions(
    chr_ref_sim_interaction_dict,
    input_output_path,
    dataset,
    method,
    data_mapping_dict,
    interaction_type,
    resolution,
    hichip_type=None,
):
    if interaction_type == "HiChIP":
        interaction_type = "HiChIP_" + hichip_type
    return np.save(
        os.path.join(
            input_output_path,
            "_".join(
                [
                    dataset,
                    method,
                    data_mapping_dict[dataset],
                    "ref_sim_average",
                    interaction_type,
                    "interaction",
                    str(resolution) + "bp.npy",
                ]
            ),
        ),
        chr_ref_sim_interaction_dict,
        allow_pickle=True,
    )


def load_average_interactions(
    input_output_path,
    dataset,
    method,
    data_mapping_dict,
    interaction_type,
    resolution,
    hichip_type=None,
):
    if interaction_type == "HiChIP":
        interaction_type = "HiChIP_" + hichip_type
    return np.load(
        os.path.join(
            input_output_path,
            "_".join(
                [
                    dataset,
                    method,
                    data_mapping_dict[dataset],
                    "ref_sim_average",
                    interaction_type,
                    "interaction",
                    str(resolution) + "bp.npy",
                ]
            ),
        ),
        allow_pickle=True,
    ).item()


def prepare_data_by_ranges(data_dict, ranges):
    ranges_dict = {}
    for length_range in ranges:
        ref_values, sim_values = [], []
        for chromosome, ref_sim_chr_dict in data_dict.items():
            for ref_cm_id, ref_sim_cm_stat_dict in ref_sim_chr_dict.items():
                if (
                    length_range[0]
                    <= ref_sim_cm_stat_dict["ref_cm_length"]
                    < length_range[1]
                ):
                    ref_values.append(ref_sim_cm_stat_dict["ref_hic_frequency"])
                    sim_values.append(ref_sim_cm_stat_dict["avg_sim_hic_frequency"])
        ranges_dict[length_range] = {"reference": ref_values, "simulated": sim_values}
    return ranges_dict


def load_interactions_per_chr(
    path_to_npy,
    interaction_type,
    ref_sim_str,
    dataset,
    method,
    data_mapping_dict,
    resolution,
):
    return np.load(
        os.path.join(
            path_to_npy,
            "_".join(
                [
                    dataset,
                    ref_sim_str,
                    method,
                    data_mapping_dict[dataset],
                    interaction_type,
                    "interactions",
                    str(resolution) + "bp.npy",
                ]
            ),
        ),
        allow_pickle=True,
    ).item()


def read_tracks_df(ref_sim_path, dataset, method, ref_sim_str, file_suffix):
    return pd.read_csv(
        os.path.join(
            ref_sim_path,
            "_".join(
                [
                    dataset,
                    ref_sim_str,
                    method,
                    file_suffix + ".tracks.bed",
                ]
            ),
        ),
        sep="\t",
        header=None,
    )


def get_cms_in_region(
    core_cm_path,
    method,
    dataset,
    chromosome,
    region_start,
    region_end,
    fully_in=True,
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
    elif method == "vcmtools":
        method_path = os.path.join(core_cm_path, method, "VCMs", vcm_window, dataset)
        tracks_path = os.path.join(
            method_path,
            dataset + "_all_VCMs_corrected_pvalue_" + str(pv_threshold) + ".tracks.bed",
        )
    elif method == "clomics":
        method_path = os.path.join(core_cm_path, method, n_peaks, bg_threshold, dataset)
        tracks_path = os.path.join(
            method_path, "_".join([dataset, "Clomics_CM.tracks.bed"])
        )
    tracks_df = pd.read_csv(tracks_path, sep="\t", header=None)
    if "chr" not in str(tracks_df.iloc[0, 0]):
        tracks_df.loc[:, "chr"] = "chr" + tracks_df.loc[:, "chr"].astype(str)
    if fully_in:
        subset_df = tracks_df.loc[
            (tracks_df.loc[:, "chr"] == chromosome)
            & (tracks_df.loc[:, "start"] >= region_start)
            & (tracks_df.loc[:, "end"] <= region_end),
            :,
        ].copy()
    else:
        subset_df = tracks_df.loc[
            (tracks_df.loc[:, "chr"] == chromosome)
            & (tracks_df.loc[:, "start"] <= region_end)
            & (tracks_df.loc[:, "end"] >= region_start),
            :,
        ].copy()
    return list(subset_df.loc[:, "cm_id"])


def get_cm_peak_dict(cm_peak_path, method, dataset):
    cm_peak_df = pd.read_csv(
        os.path.join(
            cm_peak_path,
            "_".join([dataset, method, "all_peaks.bed"]),
        ),
        sep="\t",
        header=None,
        usecols=[3, 4],
        names=["peak_id", "cm_id"],
    )
    aggregated_peaks_per_cm = (
        cm_peak_df.groupby("cm_id")["peak_id"].agg(list).to_frame().reset_index()
    )
    return dict(
        zip(aggregated_peaks_per_cm["cm_id"], aggregated_peaks_per_cm["peak_id"])
    )


def get_regions_around_cms_dict(cm_peak_path, method, dataset, frac=0.1):
    cm_peak_df = pd.read_csv(
        os.path.join(
            cm_peak_path,
            "_".join([dataset, method, "all_peaks.bed"]),
        ),
        sep="\t",
        header=None,
        usecols=[1, 2, 4],
        names=["start", "end", "cm_id"],
    )
    regions_around_cms = {}
    for cm_id, cm_peak_df in cm_peak_df.groupby("cm_id"):
        cm_start = cm_peak_df.loc[:, "start"].min()
        cm_end = cm_peak_df.loc[:, "end"].max()
        cm_length_fraction = int(round(frac * (cm_end - cm_start)))
        regions_around_cms[cm_id] = [
            cm_start - cm_length_fraction,
            cm_end + cm_length_fraction,
        ]
    return regions_around_cms


def subset_bed(bed_df, chromosome, start, end, fully_in=False):
    if fully_in:
        return bed_df.loc[
            (bed_df.loc[:, "chr"] == chromosome)
            & (bed_df.loc[:, "start"] >= start)
            & (bed_df.loc[:, "end"] <= end)
        ]
    else:
        return bed_df.loc[
            (bed_df.loc[:, "chr"] == chromosome)
            & (bed_df.loc[:, "start"] <= end)
            & (bed_df.loc[:, "end"] >= start)
        ]
