import os
import pandas as pd
import numpy as np


def get_peak_set_df_list(cm_content):
    list_for_df = []
    peak_mark_set = set()
    for cm_id, peak_str in list(zip(cm_content["cm_id"], cm_content["cm_peaks"])):
        for peak_id in peak_str.split(","):
            peak_mark_set.update([peak_id.split(":")[0]])
            list_for_df.append(
                [
                    peak_id.split(":")[-3],
                    int(peak_id.split(":")[-2]),
                    int(peak_id.split(":")[-1]),
                    peak_id,
                    cm_id,
                    "+",
                ]
            )
    return peak_mark_set, list_for_df


def read_CMs_and_store_peaks(
    input_path, file_path, method, dataset, peak_type, marks_dict, output_path
):
    cm_content = pd.read_csv(
        os.path.join(input_path, file_path),
        sep="\t",
        header=None,
    )
    cm_content.columns = ["cm_id", "cm_size", "cm_peaks"]

    if not "chr" in cm_content.iloc[0, 2]:
        for i in list(range(1, 23)):
            cm_content["cm_peaks"] = cm_content["cm_peaks"].str.replace(
                ":" + str(i) + ":", ":chr" + str(i) + ":"
            )
    list_for_df = []
    peak_mark_set = set()
    for cm_id, peak_str in list(zip(cm_content["cm_id"], cm_content["cm_peaks"])):
        for peak_id in peak_str.split(","):
            peak_mark_set.update([peak_id.split(":")[0]])
            list_for_df.append(
                [
                    peak_id.split(":")[-3],
                    int(peak_id.split(":")[-2]),
                    int(peak_id.split(":")[-1]),
                    peak_id,
                    cm_id,
                    "+",
                ]
            )
    cm_peaks_df = pd.DataFrame(list_for_df).sort_values([0, 1])
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    cm_peaks_df.to_csv(
        os.path.join(
            output_path,
            "_".join([dataset, method, peak_type, "peaks.bed"]),
        ),
        sep="\t",
        index=False,
        header=False,
    )
    if peak_type == "all":
        peaks_of_all_marks = pd.concat([marks_dict[mark] for mark in peak_mark_set])
        not_cm_peaks = peaks_of_all_marks[
            ~peaks_of_all_marks.iloc[:, 3].isin(cm_peaks_df.iloc[:, 3])
        ]
        not_cm_peaks.columns = ["chr", "start", "end", "pid", "strand"]
        not_cm_peaks = not_cm_peaks.sort_values(["chr", "start"])
        not_cm_peaks["gid"] = not_cm_peaks["pid"].copy()
        not_cm_peaks[["chr", "start", "end", "pid", "gid", "strand"]].to_csv(
            os.path.join(
                output_path,
                "_".join([dataset, "not", method, "peaks.bed"]),
            ),
            sep="\t",
            index=False,
            header=False,
        )


def get_peaks_from_totem_CMs(
    input_path, method, dataset, all_file_path, not_totem_file_path, output_path
):
    all_cm_content = pd.read_csv(
        os.path.join(input_path, all_file_path),
        sep="\t",
        header=None,
    )
    not_totem_cm_content = pd.read_csv(
        os.path.join(input_path, not_totem_file_path),
        sep="\t",
        header=None,
    )

    all_cm_content.columns = ["cm_id", "cm_size", "cm_peaks"]
    not_totem_cm_content.columns = ["cm_id", "cm_size", "cm_peaks"]

    if not "chr" in all_cm_content.iloc[0, 2]:
        for i in range(1, 23):
            all_cm_content["cm_peaks"] = all_cm_content["cm_peaks"].str.replace(
                ":" + str(i) + ":", ":chr" + str(i) + ":"
            )
    if not "chr" in not_totem_cm_content.iloc[0, 2]:
        for i in range(1, 23):
            not_totem_cm_content["cm_peaks"] = not_totem_cm_content[
                "cm_peaks"
            ].str.replace(":" + str(i) + ":", ":chr" + str(i) + ":")

    _, all_list_for_df = get_peak_set_df_list(all_cm_content)
    _, not_totem_list_for_df = get_peak_set_df_list(not_totem_cm_content)

    all_cm_peaks_df = pd.DataFrame(all_list_for_df).sort_values([0, 1])
    not_totem_cm_peaks_df = pd.DataFrame(not_totem_list_for_df).sort_values([0, 1])
    totem_cm_peaks_df = all_cm_peaks_df[
        ~all_cm_peaks_df.iloc[:, 4].isin(not_totem_cm_peaks_df.iloc[:, 4])
    ]
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    totem_cm_peaks_df.to_csv(
        os.path.join(output_path, "_".join([dataset, method, "TOTEM_peaks.bed"])),
        sep="\t",
        index=False,
        header=False,
    )


def save_npy_dict(overlap_path, overlapping_peaks_path, method, dataset, prefix=""):
    overlap_df = pd.read_csv(overlap_path, sep="\t", header=None)
    no_self_overlap_df = overlap_df[overlap_df[3] != overlap_df[9]]
    k27ac_in_third_column_df = no_self_overlap_df[
        no_self_overlap_df[3].str.startswith("H3K27ac")
    ]
    k27ac_in_third_col_peak_dict = {
        k27ac_peak_3: list(k27ac_peak_df_3[9])
        for k27ac_peak_3, k27ac_peak_df_3 in k27ac_in_third_column_df.groupby(3)
    }
    if not os.path.exists(
        os.path.join(overlapping_peaks_path, "overlapping_marks_peak_mapping", method)
    ):
        os.makedirs(
            os.path.join(
                overlapping_peaks_path, "overlapping_marks_peak_mapping", method
            )
        )
    np.save(
        os.path.join(
            overlapping_peaks_path,
            "overlapping_marks_peak_mapping",
            method,
            "_".join(
                [
                    dataset,
                    "overlapping",
                    prefix + method,
                    "peaks_k27ac_to_k4me1_mapping.npy",
                ]
            ),
        ),
        k27ac_in_third_col_peak_dict,
        allow_pickle=True,
    )
