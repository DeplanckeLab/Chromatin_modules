import os
import numpy as np
import pandas as pd


def read_bed(
    path_to_file_dir,
    file_name,
    sep="\t",
    add_chr=True,
    add_chr_col="chr",
    columns=None,
    header=None,
):
    df = pd.read_csv(os.path.join(path_to_file_dir, file_name), sep=sep, header=header)
    if columns is not None:
        df.columns = columns
    if add_chr:
        if type(add_chr_col) == str:
            if not "chr" in str(df.loc[0, add_chr_col]):
                df.loc[:, add_chr_col] = "chr" + df.loc[0, add_chr_col].astype(str)
        else:
            if not "chr" in str(df.iloc[0, add_chr_col]):
                df.iloc[:, add_chr_col] = "chr" + df.iloc[0, add_chr_col].astype(str)
    return df


def read_count_mtx(
    count_mtx_path, file_name, mark, chromosomes, add_marks=True, add_chr=True
):
    count_df = pd.read_csv(
        os.path.join(
            count_mtx_path,
            file_name,
        ),
        sep="\t",
    )
    if add_marks:
        count_df.loc[:, "mark"] = mark

    count_df = count_df.rename(columns={"#Chr": "chr"})
    if add_chr:
        # Change chromosome format if needed ('chr1', 'chr2', ...)
        if not "chr" in str(count_df.iloc[0, 0]):
            count_df.loc[:, "chr"] = "chr" + count_df.loc[:, "chr"].astype(str)

        if not "chr" in str(count_df.loc[0, "pid"]):
            count_df.loc[:, "pid"] = "chr" + count_df.loc[:, "pid"]

        if not "chr" in str(count_df.loc[0, "gid"]):
            count_df.loc[:, "gid"] = "chr" + count_df.loc[:, "gid"]

    if add_marks:
        count_df.loc[:, "pid"] = count_df.loc[:, "mark"] + ":" + count_df.loc[:, "pid"]
        count_df.drop("mark", axis=1, inplace=True)

    # Subset only specified chromosomes
    count_df = count_df.loc[count_df.loc[:, "chr"].isin(chromosomes), :].copy()
    return count_df


def get_peak_coord_in_cms_by_chr(
    cm_peaks_df, chr_col, peak_start_col, peak_end_col, cm_id_col, cm_ids_to_subset=None
):
    if cm_ids_to_subset is None:
        # 'chr1': {'cm_1': [[start_1, end_1], [start_2, end_2], ...], }
        cm_peak_coordinates_by_chr_dict = {
            chromosome: {
                cm_id: list(
                    zip(cm_df.loc[:, peak_start_col], cm_df.loc[:, peak_end_col])
                )
                for cm_id, cm_df in cm_peaks_df.loc[
                    cm_peaks_df.loc[:, chr_col] == chromosome, :
                ].groupby(cm_id_col)
            }
            for chromosome in set(cm_peaks_df.loc[:, chr_col])
        }
    else:
        cm_peak_coordinates_by_chr_dict = {
            chromosome: {
                cm_id: list(
                    zip(cm_df.loc[:, peak_start_col], cm_df.loc[:, peak_end_col])
                )
                for cm_id, cm_df in cm_peaks_df.loc[
                    cm_peaks_df.loc[:, chr_col] == chromosome, :
                ].groupby(cm_id_col)
                if cm_id in cm_ids_to_subset
            }
            for chromosome in set(cm_peaks_df.loc[:, chr_col])
        }
    return cm_peak_coordinates_by_chr_dict


def get_peak_coord_by_chr(
    peak_df,
    chr_col,
    peak_id_col,
    peak_start_col,
    peak_end_col,
):
    # 'chr1': {'H3K27ac:chr1:start:end': [start, end], ...}
    peak_coordinates_by_chr_dict = {
        chromosome: dict(
            zip(
                chr_peak_df.loc[:, peak_id_col],
                list(
                    zip(
                        chr_peak_df.loc[:, peak_start_col],
                        chr_peak_df.loc[:, peak_end_col],
                    )
                ),
            )
        )
        for chromosome, chr_peak_df in peak_df.groupby(chr_col)
    }
    return peak_coordinates_by_chr_dict


def get_peak_ids_by_cm(
    peak_df,
    peak_id_col,
    cm_id_col,
):
    # 'cm1': ['H3K27ac:chr1:start:end', ...]
    cm_peaks_dict = {
        cm_id: list(
            cm_peak_df.loc[:, peak_id_col],
        )
        for cm_id, cm_peak_df in peak_df.groupby(cm_id_col)
    }
    return cm_peaks_dict


def get_ab_compartment_dfs(df, col_comp_annotation):
    a_peaks_df = df.loc[df.loc[:, col_comp_annotation].str.contains("A")].copy()
    b_peaks_df = df.loc[df.loc[:, col_comp_annotation].str.contains("B")].copy()
    assert a_peaks_df.shape[0] + b_peaks_df.shape[0] == df.shape[0]
    return a_peaks_df, b_peaks_df


def get_cm_tracks_content_from_dict(tracks_file_path, content_file_path):
    tracks_df = pd.read_csv(tracks_file_path, sep="\t", header=None)
    content_df = pd.read_csv(content_file_path, sep="\t", header=None)
    content_df.columns = ["cm_id", "cm_size", "cm_peaks"]

    if not "chr" in content_df.loc[0, "cm_peaks"]:
        for i in range(1, 23):
            content_df.loc[:, "cm_peaks"] = content_df.loc[:, "cm_peaks"].str.replace(
                ":" + str(i) + ":", ":chr" + str(i) + ":"
            )
    return tracks_df, content_df


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


def get_dict_with_cm_tracks_df(
    path_to_filtered_cms, method, dataset, cm_type, cm_set=None, use_tad_ab_info=True
):
    if (cm_type == "ref") or (cm_type == "reference"):
        if use_tad_ab_info:
            column_names_cms = [
                "chr",
                "start",
                "end",
                "cm_id",
                "number",
                "strand",
                "start_duplicate",
                "end_duplicate",
                "numbers",
                "cm_size",
                "peak_length",
                "peak_starts",
                "is_totem",
                "chr_tad",
                "start_tad",
                "end_tad",
                "tad_id",
                "tad_cm_ovrlp_len",
                "chr_ab",
                "ab_start",
                "ab_end",
                "AB",
                "ab_cm_overlap_len",
            ]
        else:
            column_names_cms = [
                "chr",
                "start",
                "end",
                "cm_id",
                "number",
                "strand",
                "start_duplicate",
                "end_duplicate",
                "numbers",
                "cm_size",
                "peak_length",
                "peak_starts",
                "is_totem",
            ]

    else:
        column_names_cms = [
            "chr",
            "start",
            "end",
            "cm_id",
            "number",
            "strand",
            "start_duplicate",
            "end_duplicate",
            "numbers",
            "cm_size",
            "peak_length",
            "peak_starts",
            "is_totem",
            "ref_cm_id",
            "sim_cm_id",
        ]

    path_to_cms = os.path.join(
        path_to_filtered_cms,
        method,
        "all_simulated_cms",
        "_".join([dataset, cm_type, method, "no_repeating_peaks_matched"]),
    )
    tracks_df_bed = pd.read_csv(path_to_cms + ".tracks.bed", sep="\t", header=None)
    tracks_df_bed.columns = column_names_cms
    tracks_df_bed.loc[:, "cm_id"] = tracks_df_bed.loc[:, "cm_id"].str.replace("_", "~")
    if "chr" not in str(tracks_df_bed.loc[0, "chr"]):
        tracks_df_bed.loc[:, "chr"] = "chr" + tracks_df_bed.loc[:, "chr"].astype(str)
    if cm_type == "ref":
        tracks_df_bed = tracks_df_bed.loc[
            tracks_df_bed.loc[:, "cm_id"].isin(list(cm_set)), :
        ]
    dict_with_chr_cms = {
        chromosome_id: list(df_chr.loc[:, "cm_id"])
        for chromosome_id, df_chr in tracks_df_bed.groupby("chr", as_index=False)
    }
    if cm_type == "sim":
        sim_cm_set = set(tracks_df_bed["ref_cm_id"])

    if cm_type == "sim":
        return tracks_df_bed, dict_with_chr_cms, sim_cm_set
    else:
        return tracks_df_bed, dict_with_chr_cms


def get_cm_content_df(path_to_filtered_cms, method, dataset, cm_type):
    path_to_cms = os.path.join(
        path_to_filtered_cms,
        method,
        "all_simulated_cms",
        "_".join([dataset, cm_type, method, "no_repeating_peaks_matched"]),
    )
    content_df_bed = pd.read_csv(path_to_cms + ".content.txt", sep="\t", header=None)
    content_df_bed.columns = ["cm_id", "n_peaks", "peaks"]
    content_df_bed.loc[:, "cm_id"] = content_df_bed.loc[:, "cm_id"].str.replace(
        "_", "~"
    )
    return content_df_bed


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


def get_merged_ref_sim(dict_with_ref_sim_paths, simulation_type, use_tad_ab_info=True):
    reference_file = dict_with_ref_sim_paths[simulation_type]["reference"]
    simulated_file = dict_with_ref_sim_paths[simulation_type]["simulated"]

    reference_tracks = pd.read_csv(
        reference_file + ".tracks.bed",
        sep="\t",
        header=None,
    )
    simulated_tracks = pd.read_csv(
        simulated_file + ".tracks.bed",
        sep="\t",
        header=None,
    )
    if use_tad_ab_info:
        ref_cols = [
            "chr",
            "start",
            "end",
            "cm_id",
            "number",
            "strand",
            "start_duplicate",
            "end_duplicate",
            "numbers",
            "cm_size",
            "peak_lengths",
            "peak_starts",
            "is_totem",
            "chr_tad",
            "start_tad",
            "end_tad",
            "tad_id",
            "tad_cm_ovrlp_len",
            "chr_ab",
            "ab_start",
            "ab_end",
            "AB",
            "ab_cm_overlap_len",
        ]
    else:
        ref_cols = [
            "chr",
            "start",
            "end",
            "cm_id",
            "number",
            "strand",
            "start_duplicate",
            "end_duplicate",
            "numbers",
            "cm_size",
            "peak_lengths",
            "peak_starts",
            "is_totem",
        ]

    sim_cols = [
        "chr",
        "start",
        "end",
        "cm_id",
        "number",
        "strand",
        "start_duplicate",
        "end_duplicate",
        "numbers",
        "cm_size",
        "peak_lengths",
        "peak_starts",
        "is_totem",
    ]
    if simulation_type == "all_simulated_cms":
        reference_tracks.columns = ["ref_" + r_col_name for r_col_name in ref_cols]
        simulated_tracks.columns = ["sim_" + s_col_name for s_col_name in sim_cols] + [
            "ref_cm_id",
            "sim_ref_cm_id",
        ]
    else:
        reference_tracks.columns = ["ref_" + r_col_name for r_col_name in ref_cols] + [
            "ref_some_num"
        ]
        simulated_tracks.columns = (
            ["sim_" + s_col_name for s_col_name in sim_cols]
            + ["ref_cm_id", "sim_ref_cm_id"]
            + ["sim_some_num"]
        )

    reference_tracks.loc[:, "ref_cm_id"] = reference_tracks.loc[
        :, "ref_cm_id"
    ].str.replace("_", "~")
    if not "chr" in str(reference_tracks.iloc[0, 0]):
        reference_tracks.loc[:, "ref_chr"] = "chr" + reference_tracks.loc[
            :, "ref_chr"
        ].astype(str)

    if not "chr" in str(simulated_tracks.iloc[0, 0]):
        simulated_tracks.loc[:, "sim_chr"] = "chr" + simulated_tracks.loc[
            :, "sim_chr"
        ].astype(str)

    reference_tracks.loc[:, "ref_length"] = (
        reference_tracks.loc[:, "ref_end"] - reference_tracks.loc[:, "ref_start"]
    )
    simulated_tracks.loc[:, "sim_length"] = (
        simulated_tracks.loc[:, "sim_end"] - simulated_tracks.loc[:, "sim_start"]
    )
    reference_tracks.loc[:, "ref_peak_bp_length"] = [
        np.sum([int(ref_peak_length) for ref_peak_length in ref_p_lengths])
        for ref_p_lengths in reference_tracks.loc[:, "ref_peak_lengths"]
        .str.split(",")
        .tolist()
    ]
    simulated_tracks.loc[:, "sim_peak_bp_length"] = [
        np.sum([int(sim_peak_length) for sim_peak_length in sim_p_lengths])
        for sim_p_lengths in simulated_tracks.loc[:, "sim_peak_lengths"]
        .str.split(",")
        .tolist()
    ]
    ref_sim_merged = simulated_tracks.merge(reference_tracks, on="ref_cm_id").loc[
        :,
        [
            "ref_chr",
            "sim_chr",
            "ref_cm_id",
            "sim_cm_id",
            "ref_length",
            "sim_length",
            "ref_peak_bp_length",
            "sim_peak_bp_length",
        ],
    ]
    return ref_sim_merged


def get_gc_content_per_cm(cm_peaks_dict, gc_content_dict):
    return {
        cm_id: np.mean([gc_content_dict.get(cm_peak) for cm_peak in cm_peaks])
        for cm_id, cm_peaks in cm_peaks_dict.items()
    }


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


def check_chr(data_frame, col_idx, data_type):
    if not "chr" in data_frame.iloc[0, col_idx]:
        if data_type == "content":
            for i in range(1, 23):
                data_frame.iloc[:, col_idx] = data_frame.iloc[:, col_idx].str.replace(
                    ":" + str(i) + ":", ":chr" + str(i) + ":"
                )
        elif data_type == "tracks":
            data_frame.iloc[:, col_idx] = (
                data_frame.iloc[:, col_idx].astype(str).str.replace("chr", "")
            )
            data_frame.iloc[:, col_idx] = "chr" + data_frame.iloc[:, col_idx].astype(
                str
            )
    return data_frame
