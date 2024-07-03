import os
import collections
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib_venn
import upsetplot


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
        method_path = os.path.join(
            core_cm_path, method, "vcmtoolss", vcm_window, dataset
        )
        tracks_path = os.path.join(
            method_path,
            dataset
            + "_all_vcmtoolss_corrected_pvalue_"
            + str(pv_threshold)
            + ".tracks.bed",
        )
        content_path = os.path.join(
            method_path,
            dataset
            + "_all_vcmtoolss_corrected_pvalue_"
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


def get_similarity_scores_as_multiset_df(ds_pair_scores_dict, dataset, ds_type="first"):
    if ds_type == "first":
        cm_idx = 0
    else:
        cm_idx = 1
    ds_similarity_list = [
        cm_pair_score
        for cm_pair, cm_pair_score in ds_pair_scores_dict.items()
        if cm_pair.split("_")[cm_idx] != "NAN"
    ]
    multiplicities_dict = dict(collections.Counter(ds_similarity_list))
    assert sum(multiplicities_dict.values()) == len(ds_similarity_list)
    return pd.DataFrame.from_dict(
        multiplicities_dict, orient="index", columns=[dataset]
    )


def process_thresholds(cm_overlap, method, dataset, path, threshold, min_cm_size):
    if threshold == 0:
        return cm_overlap.get_unique_CMs_df(
            method,
            dataset,
            path,
            threshold,
            min_cm_size=min_cm_size,
            save_filtered_tracks=True,
        )
    else:
        return cm_overlap.get_universal_CMs_df(
            method,
            dataset,
            path,
            threshold,
            min_cm_size=min_cm_size,
            save_filtered_tracks=True,
            direction="both",
        )


def handle_threshold(threshold_str):
    if "-" in threshold_str:
        range_start, rest = threshold_str.split("-")
        range_end, step = rest.split(";")
        return np.arange(
            float(range_start), float(range_end) + float(step), float(step)
        )
    else:
        return [float(threshold_str)]


def get_sets_plot_venn(
    cm_pairs_for_methods, dataset, total_pairs, threshold, output_path
):
    method_1 = "clomics"
    method_2 = "vcmtools"
    method_3 = "phm"
    unique_crds_1 = cm_pairs_for_methods[method_1 + "_" + method_2]["mthd_1_unique"]
    unique_vcms_1 = cm_pairs_for_methods[method_1 + "_" + method_2]["mthd_2_unique"]

    unique_crds_2 = cm_pairs_for_methods[method_1 + "_" + method_3]["mthd_1_unique"]
    unique_phms_1 = cm_pairs_for_methods[method_1 + "_" + method_3]["mthd_2_unique"]

    unique_vcms_2 = cm_pairs_for_methods[method_2 + "_" + method_3]["mthd_1_unique"]
    unique_phms_2 = cm_pairs_for_methods[method_2 + "_" + method_3]["mthd_2_unique"]

    unique_crds = unique_crds_1.intersection(unique_crds_2)
    unique_vcms = unique_vcms_1.intersection(unique_vcms_2)
    unique_phms = unique_phms_1.intersection(unique_phms_2)

    ovrlp_crds_1 = cm_pairs_for_methods[method_1 + "_" + method_2][
        "mthd_1_specific_cms"
    ]
    ovrlp_vcms_1 = cm_pairs_for_methods[method_1 + "_" + method_2][
        "mthd_2_specific_cms"
    ]

    ovrlp_crds_2 = cm_pairs_for_methods[method_1 + "_" + method_3][
        "mthd_1_specific_cms"
    ]
    ovrlp_phms_1 = cm_pairs_for_methods[method_1 + "_" + method_3][
        "mthd_2_specific_cms"
    ]

    ovrlp_vcms_2 = cm_pairs_for_methods[method_2 + "_" + method_3][
        "mthd_1_specific_cms"
    ]
    ovrlp_phms_2 = cm_pairs_for_methods[method_2 + "_" + method_3][
        "mthd_2_specific_cms"
    ]

    common_crds = ovrlp_crds_1.intersection(ovrlp_crds_2)
    common_vcms = ovrlp_vcms_1.intersection(ovrlp_vcms_2)
    common_phms = ovrlp_phms_1.intersection(ovrlp_phms_2)
    common_cms = {"clomics": common_crds, "vcmtools": common_vcms, "phm": common_phms}
    method_pairs = [["clomics", "vcmtools"], ["clomics", "phm"], ["vcmtools", "phm"]]
    # Common CMs in all methods:
    cms_common_in_all_methods = set()
    for method_1, method_2 in method_pairs:
        for cm_pair in cm_pairs_for_methods[method_1 + "_" + method_2]["overlapping"]:
            cm_1, cm_2 = cm_pair.split("_")
            if (cm_1 in common_cms[method_1]) and (cm_2 in common_cms[method_2]):
                cms_common_in_all_methods.update([cm_pair])

    # CMs common in method pairs
    cm_sets_for_method_pairs = {}
    for method_1, method_2 in method_pairs:
        overlapping_cms = cm_pairs_for_methods[method_1 + "_" + method_2]["overlapping"]
        method_1_method_2_only_cms = overlapping_cms - cms_common_in_all_methods
        cm_sets_for_method_pairs[method_1 + "_" + method_2] = method_1_method_2_only_cms

    crd_vcm_common = cm_sets_for_method_pairs["clomics_vcmtools"]
    crd_phm_common = cm_sets_for_method_pairs["clomics_phm"]
    vcm_phm_common = cm_sets_for_method_pairs["vcmtools_phm"]
    vcm_in_pairs = {
        cm for pair in vcm_phm_common for cm in pair.split("_") if cm.startswith("vcm")
    }.union(
        {
            cm
            for pair in crd_vcm_common
            for cm in pair.split("_")
            if cm.startswith("vcm")
        }
    )
    phm_in_pairs = {
        cm for pair in vcm_phm_common for cm in pair.split("_") if cm.startswith("phm")
    }.union(
        {
            cm
            for pair in crd_phm_common
            for cm in pair.split("_")
            if cm.startswith("phm")
        }
    )
    crd_in_pairs = {
        cm for pair in crd_phm_common for cm in pair.split("_") if cm.startswith("crd")
    }.union(
        {
            cm
            for pair in crd_vcm_common
            for cm in pair.split("_")
            if cm.startswith("crd")
        }
    )

    common_cms_all = {
        cm for pair in cms_common_in_all_methods for cm in pair.split("_")
    }
    unique_cms = {
        "clomics": unique_crds - common_cms_all - crd_in_pairs,
        "vcmtools": unique_vcms - common_cms_all - vcm_in_pairs,
        "phm": unique_phms - common_cms_all - phm_in_pairs,
    }
    total_set_size = (
        len(cms_common_in_all_methods)
        + len(unique_cms["clomics"])
        + len(unique_cms["vcmtools"])
        + len(unique_cms["phm"])
        + len(crd_vcm_common)
        + len(crd_phm_common)
        + len(vcm_phm_common)
    )
    total_n_cm_pairs = (
        total_pairs
        + len(unique_cms["clomics"])
        + len(unique_cms["vcmtools"])
        + len(unique_cms["phm"])
    )
    assert total_n_cm_pairs == total_set_size
    plt.figure(figsize=(5, 5))
    venn_d = matplotlib_venn.venn3(
        subsets=(
            len(unique_cms["vcmtools"]),
            len(unique_cms["clomics"]),
            len(crd_vcm_common),
            len(unique_cms["phm"]),
            len(vcm_phm_common),
            len(crd_phm_common),
            len(cms_common_in_all_methods),
        ),
        set_labels=("vcmtools", "clomics", "phm"),
        set_colors=("#e8c6b0", "#b0b0b2", "#80b0ae"),
        subset_label_formatter=lambda x: f"{x/total_n_cm_pairs:.1%}",
        alpha=0.7,
    )

    # outline of circle line style and width
    matplotlib_venn.venn3_circles(
        subsets=(
            len(unique_cms["vcmtools"]),
            len(unique_cms["clomics"]),
            len(crd_vcm_common),
            len(unique_cms["phm"]),
            len(vcm_phm_common),
            len(crd_phm_common),
            len(cms_common_in_all_methods),
        ),
        #     linestyle='dashed',
        #     ls='--',
        linewidth=0.5,
    )

    for text in venn_d.set_labels:
        text.set_fontsize(12)

    plt.title(
        "\n".join(
            [
                "Methods overlap for CM pairs in " + dataset + ",",
                "total #CM pairs=" + str(total_n_cm_pairs) + ",",
                "similarity threshold for uniqueness = " + str(threshold) + "\n",
            ]
        ),
        size=12,
    )

    plt.savefig(
        os.path.join(
            output_path,
            "_".join(
                [
                    "methods_overlap_as_CM_pairs",
                    dataset,
                    "similarity_threshold_for_uniqueness",
                    str(threshold) + ".pdf",
                ]
            ),
        ),
        dpi=300,
        transparent=True,
        bbox_inches="tight",
    )


def get_cm_membership_df(cm_pairs_for_methods, total_pairs):
    method_1 = "clomics"
    method_2 = "vcmtools"
    method_3 = "phm"
    unique_crds_1 = cm_pairs_for_methods[method_1 + "_" + method_2]["mthd_1_unique"]
    unique_vcms_1 = cm_pairs_for_methods[method_1 + "_" + method_2]["mthd_2_unique"]

    unique_crds_2 = cm_pairs_for_methods[method_1 + "_" + method_3]["mthd_1_unique"]
    unique_phms_1 = cm_pairs_for_methods[method_1 + "_" + method_3]["mthd_2_unique"]

    unique_vcms_2 = cm_pairs_for_methods[method_2 + "_" + method_3]["mthd_1_unique"]
    unique_phms_2 = cm_pairs_for_methods[method_2 + "_" + method_3]["mthd_2_unique"]

    unique_crds = unique_crds_1.intersection(unique_crds_2)
    unique_vcms = unique_vcms_1.intersection(unique_vcms_2)
    unique_phms = unique_phms_1.intersection(unique_phms_2)

    ovrlp_crds_1 = cm_pairs_for_methods[method_1 + "_" + method_2][
        "mthd_1_specific_cms"
    ]
    ovrlp_vcms_1 = cm_pairs_for_methods[method_1 + "_" + method_2][
        "mthd_2_specific_cms"
    ]

    ovrlp_crds_2 = cm_pairs_for_methods[method_1 + "_" + method_3][
        "mthd_1_specific_cms"
    ]
    ovrlp_phms_1 = cm_pairs_for_methods[method_1 + "_" + method_3][
        "mthd_2_specific_cms"
    ]

    ovrlp_vcms_2 = cm_pairs_for_methods[method_2 + "_" + method_3][
        "mthd_1_specific_cms"
    ]
    ovrlp_phms_2 = cm_pairs_for_methods[method_2 + "_" + method_3][
        "mthd_2_specific_cms"
    ]

    common_crds = ovrlp_crds_1.intersection(ovrlp_crds_2)
    common_vcms = ovrlp_vcms_1.intersection(ovrlp_vcms_2)
    common_phms = ovrlp_phms_1.intersection(ovrlp_phms_2)
    common_cms = {"clomics": common_crds, "vcmtools": common_vcms, "phm": common_phms}
    method_pairs = [["clomics", "vcmtools"], ["clomics", "phm"], ["vcmtools", "phm"]]
    # Common CMs in all methods:
    cms_common_in_all_methods = set()
    for method_1, method_2 in method_pairs:
        for cm_pair in cm_pairs_for_methods[method_1 + "_" + method_2]["overlapping"]:
            cm_1, cm_2 = cm_pair.split("_")
            if (cm_1 in common_cms[method_1]) and (cm_2 in common_cms[method_2]):
                cms_common_in_all_methods.update([cm_pair])

    # CMs common in method pairs
    cm_sets_for_method_pairs = {}
    for method_1, method_2 in method_pairs:
        overlapping_cms = cm_pairs_for_methods[method_1 + "_" + method_2]["overlapping"]
        method_1_method_2_only_cms = overlapping_cms - cms_common_in_all_methods
        cm_sets_for_method_pairs[method_1 + "_" + method_2] = method_1_method_2_only_cms

    crd_vcm_common = cm_sets_for_method_pairs["clomics_vcmtools"]
    crd_phm_common = cm_sets_for_method_pairs["clomics_phm"]
    vcm_phm_common = cm_sets_for_method_pairs["vcmtools_phm"]
    vcm_in_pairs = {
        cm for pair in vcm_phm_common for cm in pair.split("_") if cm.startswith("vcm")
    }.union(
        {
            cm
            for pair in crd_vcm_common
            for cm in pair.split("_")
            if cm.startswith("vcm")
        }
    )
    phm_in_pairs = {
        cm for pair in vcm_phm_common for cm in pair.split("_") if cm.startswith("phm")
    }.union(
        {
            cm
            for pair in crd_phm_common
            for cm in pair.split("_")
            if cm.startswith("phm")
        }
    )
    crd_in_pairs = {
        cm for pair in crd_phm_common for cm in pair.split("_") if cm.startswith("crd")
    }.union(
        {
            cm
            for pair in crd_vcm_common
            for cm in pair.split("_")
            if cm.startswith("crd")
        }
    )

    common_cms_all = {
        cm for pair in cms_common_in_all_methods for cm in pair.split("_")
    }
    unique_cms = {
        "clomics": unique_crds - common_cms_all - crd_in_pairs,
        "vcmtools": unique_vcms - common_cms_all - vcm_in_pairs,
        "phm": unique_phms - common_cms_all - phm_in_pairs,
    }
    total_set_size = (
        len(cms_common_in_all_methods)
        + len(unique_cms["clomics"])
        + len(unique_cms["vcmtools"])
        + len(unique_cms["phm"])
        + len(crd_vcm_common)
        + len(crd_phm_common)
        + len(vcm_phm_common)
    )
    total_n_cm_pairs = (
        total_pairs
        + len(unique_cms["clomics"])
        + len(unique_cms["vcmtools"])
        + len(unique_cms["phm"])
    )
    assert total_n_cm_pairs == total_set_size
    # print("total n CM pairs", total_n_cm_pairs)

    cms_group_assignments = {
        "vcmtools": unique_cms["vcmtools"],
        "clomics": unique_cms["clomics"],
        "phm": unique_cms["phm"],
        "clomics~vcmtools": crd_vcm_common,
        "vcmtools~phm": vcm_phm_common,
        "clomics~phm": crd_phm_common,
        "clomics~vcmtools~phm": cms_common_in_all_methods,
    }
    cm_memership_list = [
        group.split("~")
        for group, cms_in_group in cms_group_assignments.items()
        for cm in cms_in_group
    ]

    return upsetplot.from_memberships(cm_memership_list)
