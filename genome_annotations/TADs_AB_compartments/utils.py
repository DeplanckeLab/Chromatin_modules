import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scipy

matplotlib.rcParams["axes.spines.right"] = False
matplotlib.rcParams["axes.spines.top"] = False


def boxplot_by_range(
    dict_with_sth_per_method,
    characteristic,
    dataset,
    save_figure=False,
    path=None,
    file_name=None,
    stars_y_shift=None,
):
    color_dict = {"vcmtools": "#E8C6B0", "clomics": "#B0B0B2", "phm": "#80B0AE"}
    labels = list(dict_with_sth_per_method.keys())
    x = np.arange(1, len(labels) * 2 + 1, 2)
    positions = list(np.ravel(list(zip(x, x + 0.5))))
    x_ticks = x + (0.5 / 2)

    values_by_dataset, significant_differences, stars = get_values_by_dataset(
        dict_with_sth_per_method, x, x + 0.5
    )
    fig, ax = plt.subplots(1, 1, figsize=(3, 5))
    bplot = ax.violinplot(
        values_by_dataset, positions=positions, widths=0.4, showmedians=True
    )
    trans = matplotlib.transforms.blended_transform_factory(ax.transData, ax.transData)
    for s, (start, end) in enumerate(significant_differences):
        left_array = values_by_dataset[positions.index(start)]
        right_array = values_by_dataset[positions.index(end)]

        max_1 = np.quantile(left_array, 0.75) + 1.5 * (
            np.quantile(left_array, 0.75) - np.quantile(left_array, 0.25)
        )
        sorted_lst_1 = np.array(sorted(left_array))[sorted(left_array) <= max_1]
        whisker_1 = sorted_lst_1[len(sorted_lst_1) - 1]

        max_2 = np.quantile(right_array, 0.75) + 1.5 * (
            np.quantile(right_array, 0.75) - np.quantile(right_array, 0.25)
        )

        sorted_lst_2 = np.array(sorted(right_array))[sorted(right_array) <= max_2]
        whisker_2 = sorted_lst_2[len(sorted_lst_2) - 1]
        if characteristic == "length":
            q75 = max(whisker_1, whisker_2) + stars_y_shift
        else:
            q75 = 7.5
        ax.annotate(
            stars[s],
            xy=(start + 0.25, q75),
            xytext=(start + 0.25, q75),
            transform=trans,
            ha="center",
            va="bottom",
            arrowprops=dict(
                arrowstyle="-[, widthB="
                + str((end - start) * 2 - 0.05)
                + ", lengthB=0.5",
                lw=0.7,
            ),
        )

    for partname in ("cbars", "cmedians", "cmins", "cmaxes"):
        bplot[partname].set_edgecolor("black")
        bplot[partname].set_linewidth(1)
    labels_for_colors = np.ravel(list(zip(labels, labels)))
    for j, box in enumerate(bplot["bodies"]):
        box.set(
            facecolor=color_dict[labels_for_colors[j]],
            edgecolor="black",
            linewidth=1,
        )
        if j % 2 == 0:
            box.set_alpha(0.5)
    ax.set_title("Chromatin module " + characteristic + " w.r.t. A/B compartments\n")

    label_patches = [
        matplotlib.patches.Patch(color="gray", alpha=1, label="A compartment"),
        matplotlib.patches.Patch(color="lightgray", alpha=0.5, label="B compartment"),
    ]

    # Put a legend below current axis
    ax.legend(
        handles=label_patches,
        loc="upper center",
        bbox_to_anchor=(1.7, 0.5),
        fancybox=False,
        shadow=False,
        ncol=1,
    )
    label_names = []
    for label in labels:
        if label == "vcmtools":
            label_names.append("VCMtools")
        if label == "clomics":
            label_names.append("Clomics")
        elif label == "phm":
            label_names.append("PHM")
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(label_names, rotation=0, size=12)

    if characteristic == "length":
        ax.set_yticks(np.arange(2, 9))
        ax.set_yticklabels(["$10^" + str(i) + "$" for i in np.arange(2, 9)])
        ax.set_ylim((2, 8))
        ax.set_ylabel("Chromatin module " + characteristic + " in bp", size=12)
    else:
        ax.set_yticks(np.arange(1, 8))
        ax.set_yticklabels([2**i for i in np.arange(1, 8)])
        ax.set_ylim((1, 8))
        ax.set_ylabel("Chromatin module " + characteristic, size=12)
    ax.set_title(dataset)
    if save_figure:
        if not os.path.exists(path):
            os.makedirs(path)
        plt.savefig(
            os.path.join(path, file_name),
            dpi=300,
            transparent=True,
            bbox_inches="tight",
        )


def get_values_by_dataset(dataset_dict, a, b):
    values_by_dataset = []
    significant_differences = []
    stars = []
    for pos, (method, AB_dict) in enumerate(dataset_dict.items()):
        lst_1 = np.array(AB_dict["A"])
        lst_2 = np.array(AB_dict["B"])
        values_by_dataset.append(lst_1)
        values_by_dataset.append(lst_2)

        if (len(lst_1) < 2) or (len(lst_2) < 2):
            pv = 1
        else:
            _, pv = scipy.stats.mannwhitneyu(lst_1, lst_2)
        if pv <= 0.05:
            significant_differences.append([a[pos], b[pos]])
            if pv <= 0.05:
                star = "*"
            if pv <= 0.01:
                star = "**"
            if pv <= 0.001:
                star = "***"
            if pv <= 0.0001:
                star = "****"
            stars.append(star)
    return values_by_dataset, significant_differences, stars


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
    content_df.columns = ["cm_id", "cm_size", "cm_peaks"]

    if not "chr" in content_df.loc[0, "cm_peaks"]:
        for i in range(1, 23):
            content_df.loc[:, "cm_peaks"] = content_df.loc[:, "cm_peaks"].str.replace(
                ":" + str(i) + ":", ":chr" + str(i) + ":"
            )
    return tracks_df, content_df
