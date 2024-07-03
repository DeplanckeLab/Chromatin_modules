import gseapy as gp
import itertools
import networkx as nx
import numpy as np
import os
import pandas as pd
import tqdm.auto as tqdm
import sys
import warnings

warnings.filterwarnings("ignore")


class CM_Overlap:
    def __init__(
        self, output_path=None, methods_list=None, score_type=None, use_cmQTLs=False
    ):
        self.output_path = output_path
        self.methods_list = methods_list
        self.use_cmQTLs = use_cmQTLs
        self.score_type = score_type
        self.n_datasets = len(methods_list)
        if self.n_datasets < 2:
            sys.exit("ERROR: less than two datasets were provided. Check the input")

    @staticmethod
    def subset_cms_by_size(cm_tracks_path, min_cm_size):
        cm_df = pd.read_csv(cm_tracks_path, sep="\t", header=None)
        cm_df.iloc[:, 3] = cm_df.iloc[:, 3].str.replace("_", "~")
        return list(cm_df[cm_df.iloc[:, 9] >= min_cm_size].iloc[:, 3])

    @staticmethod
    def get_cms_with_cmQTL(core_cm_path, method, dataset):
        crd_qtls = pd.read_csv(
            os.path.join(
                core_cm_path,
                "aCM_scores",
                method,
                "_".join([method, dataset, "cmQTLs.significant.txt"]),
            ),
            sep=" ",
            header=None,
        )
        return list(crd_qtls.iloc[:, 0].str.replace("_", "~"))

    @staticmethod
    def get_cm_peaks(input_peak_path):
        peak_df = pd.read_csv(input_peak_path, sep="\t", header=None)
        peak_df.iloc[:, 4] = peak_df.iloc[:, 4].str.replace("_", "~")
        return peak_df

    @staticmethod
    def merge_intersecting_intervals(list_of_pairs):
        merged_intervals = []
        for pair in sorted(list_of_pairs):
            if not merged_intervals:
                merged_intervals.append((pair[0], pair[1]))
            else:
                lower_pair = merged_intervals[-1]
                if pair[0] <= lower_pair[1]:
                    upper_bound = max(lower_pair[1], pair[1])
                    merged_intervals[-1] = (
                        lower_pair[0],
                        upper_bound,
                    )  # replace by merged interval
                else:
                    merged_intervals.append((pair[0], pair[1]))
        return merged_intervals

    def save_tracks(
        self,
        tracks_df,
        dataset,
        method,
        similarity_threshold,
        tracks_output_path,
        direction=None,
    ):
        if (similarity_threshold > 0) and (similarity_threshold < 1):
            if direction == "greater":
                tracks_df.to_csv(
                    os.path.join(
                        tracks_output_path,
                        "_".join(
                            [
                                dataset,
                                method,
                                self.score_type,
                                "scores",
                                "greater",
                                str(similarity_threshold) + ".tracks.bed",
                            ]
                        ),
                    ),
                    sep="\t",
                    index=False,
                    header=False,
                )
            else:
                tracks_df.to_csv(
                    os.path.join(
                        tracks_output_path,
                        "_".join(
                            [
                                dataset,
                                method,
                                self.score_type,
                                "scores",
                                "less",
                                str(similarity_threshold) + ".tracks.bed",
                            ]
                        ),
                    ),
                    sep="\t",
                    index=False,
                    header=False,
                )
        elif (similarity_threshold == 0) or (similarity_threshold == 1):
            suffix = "equal"
            tracks_df.to_csv(
                os.path.join(
                    tracks_output_path,
                    "_".join(
                        [
                            dataset,
                            method,
                            self.score_type,
                            "scores",
                            suffix,
                            str(similarity_threshold) + ".tracks.bed",
                        ]
                    ),
                ),
                sep="\t",
                index=False,
                header=False,
            )
        elif similarity_threshold is None:
            suffix = "all"
            tracks_df.to_csv(
                os.path.join(
                    tracks_output_path,
                    "_".join(
                        [
                            dataset,
                            method,
                            self.score_type,
                            "scores",
                            suffix,
                            str(similarity_threshold) + ".tracks.bed",
                        ]
                    ),
                ),
                sep="\t",
                index=False,
                header=False,
            )

    def load_tracks(
        self, dataset, method, similarity_threshold, tracks_output_path, direction=None
    ):
        if (similarity_threshold > 0) and (similarity_threshold < 1):
            suffix = direction
            if suffix == "both":
                df_1 = pd.read_csv(
                    os.path.join(
                        tracks_output_path,
                        "_".join(
                            [
                                dataset,
                                method,
                                self.score_type,
                                "scores",
                                "less",
                                str(similarity_threshold) + ".tracks.bed",
                            ]
                        ),
                    ),
                    sep="\t",
                    header=None,
                )
                df_2 = pd.read_csv(
                    os.path.join(
                        tracks_output_path,
                        "_".join(
                            [
                                dataset,
                                method,
                                self.score_type,
                                "scores",
                                "greater",
                                str(similarity_threshold) + ".tracks.bed",
                            ]
                        ),
                    ),
                    sep="\t",
                    header=None,
                )
                return df_1, df_2
        elif (similarity_threshold == 0) or (similarity_threshold == 1):
            suffix = "equal"
        elif similarity_threshold is None:
            suffix = "all"
        return pd.read_csv(
            os.path.join(
                tracks_output_path,
                "_".join(
                    [
                        dataset,
                        method,
                        self.score_type,
                        "scores",
                        suffix,
                        str(similarity_threshold) + ".tracks.bed",
                    ]
                ),
            ),
            sep="\t",
            header=None,
        )

    def save_enrichr_df(
        self,
        method,
        dataset,
        similarity_threshold,
        enrichr_df,
        enrichr_library,
        pv_filtered_output,
        enrichment_pv_threshold,
        how=None,
        direction="",
    ):
        if self.use_cmQTLs:
            folder_name = "modules_with_QTL"
        else:
            folder_name = "all_modules"
        if how == "union":
            enr_output_path = os.path.join(
                self.output_path,
                "grouped_enrichment_results",
                folder_name,
                str(similarity_threshold),
            )
        else:
            enr_output_path = os.path.join(
                self.output_path,
                "enrichr_output",
                folder_name,
                str(similarity_threshold),
            )

        if not os.path.exists(enr_output_path):
            os.makedirs(enr_output_path)
        if pv_filtered_output:
            suffix = str(enrichment_pv_threshold) + "_pv_cutoff"
        else:
            suffix = "full"
        enrichr_df.to_csv(
            os.path.join(
                enr_output_path,
                "_".join(
                    [
                        method,
                        dataset,
                        direction,
                        str(similarity_threshold),
                        enrichr_library,
                        "enrichr",
                        suffix,
                        "df.bed",
                    ]
                ),
            ),
            sep="\t",
            index=False,
            header=True,
        )

    def save_enriched_genes(
        self,
        method,
        dataset,
        similarity_threshold,
        gene_df,
        enrichr_library,
        pv_filtered_output,
        enrichment_pv_threshold,
        direction="",
    ):
        if self.use_cmQTLs:
            folder_name = "modules_with_QTL"
        else:
            folder_name = "all_modules"
        enr_output_path = os.path.join(
            self.output_path, "enrichr_output", folder_name, str(similarity_threshold)
        )
        if not os.path.exists(enr_output_path):
            os.makedirs(enr_output_path)
        if pv_filtered_output:
            suffix = str(enrichment_pv_threshold) + "_pv_cutoff"
        else:
            suffix = "full"
        gene_df.to_csv(
            os.path.join(
                enr_output_path,
                "_".join(
                    [
                        method,
                        dataset,
                        direction,
                        str(similarity_threshold),
                        enrichr_library,
                        "enrichr",
                        suffix,
                        "gene_counts.txt",
                    ]
                ),
            ),
            sep="\t",
            index=False,
            header=True,
        )

    @staticmethod
    def get_available_enrichr_libraries():
        return gp.get_library_name()

    def save_npy(
        self,
        method_A,
        method_B,
        dataset,
        F1_dict,
        overlap_dict_path,
    ):
        prefix = "_".join([method_A, method_B])
        if self.use_cmQTLs:
            suffix = "scores_dict_only_with_QTL.npy"
        else:
            suffix = "scores_dict.npy"
        np.save(
            os.path.join(
                overlap_dict_path, "_".join([prefix, dataset, self.score_type, suffix])
            ),
            F1_dict,
        )

    def load_npy(self, method_A, method_B, dataset, overlap_dict_path):
        prefix = "_".join([method_A, method_B])
        if self.use_cmQTLs:
            suffix = "scores_dict_only_with_QTL.npy"
        else:
            suffix = "scores_dict.npy"
        return np.load(
            os.path.join(
                overlap_dict_path, "_".join([prefix, dataset, self.score_type, suffix])
            ),
            allow_pickle=True,
        ).item()

    def get_jaccard_peak_overlap(
        self,
        method_AB_peak_overlap_df,
        method_A_peak_start_col=None,
        method_A_peak_end_col=None,
        method_A_vcm_id_col=None,
        method_B_peak_start_col=None,
        method_B_peak_end_col=None,
        method_B_vcm_id_col=None,
        peak_len_overlap_col=None,
        A_CMs_with_qtl=None,
        B_CMs_with_qtl=None,
        bp_jaccard=True,
    ):
        if self.use_cmQTLs:
            if (A_CMs_with_qtl is None) and (B_CMs_with_qtl is None):
                sys.exit(
                    "No chromatin modules with QTLs were provided, check the input!\n Exiting..."
                )
            if (A_CMs_with_qtl is None) and (B_CMs_with_qtl is not None):
                method_AB_peak_overlap_df = method_AB_peak_overlap_df[
                    method_AB_peak_overlap_df.iloc[:, method_B_vcm_id_col].isin(
                        B_CMs_with_qtl
                    )
                ].copy()
            elif (A_CMs_with_qtl is not None) and (B_CMs_with_qtl is None):
                method_AB_peak_overlap_df = method_AB_peak_overlap_df[
                    method_AB_peak_overlap_df.iloc[:, method_A_vcm_id_col].isin(
                        A_CMs_with_qtl
                    )
                ].copy()
            else:
                method_AB_peak_overlap_df = method_AB_peak_overlap_df[
                    method_AB_peak_overlap_df.iloc[:, method_A_vcm_id_col].isin(
                        A_CMs_with_qtl
                    )
                    & method_AB_peak_overlap_df.iloc[:, method_B_vcm_id_col].isin(
                        B_CMs_with_qtl
                    )
                ].copy()

        method_AB_peak_overlap_df.loc[:, "method_A_peak_length"] = (
            method_AB_peak_overlap_df[method_A_peak_end_col].copy()
            - method_AB_peak_overlap_df[method_A_peak_start_col].copy()
        )
        method_AB_peak_overlap_df.loc[:, "method_B_peak_length"] = (
            method_AB_peak_overlap_df[method_B_peak_end_col].copy()
            - method_AB_peak_overlap_df[method_B_peak_start_col].copy()
        )
        if bp_jaccard:
            method_AB_peak_overlap_df.loc[:, "jaccard"] = method_AB_peak_overlap_df[
                peak_len_overlap_col
            ].copy() / (
                method_AB_peak_overlap_df.loc[:, "method_A_peak_length"].copy()
                + method_AB_peak_overlap_df.loc[:, "method_B_peak_length"].copy()
                - method_AB_peak_overlap_df[peak_len_overlap_col].copy()
            )
        else:
            method_AB_peak_overlap_df.loc[:, "jaccard"] = 1

        method_AB_overlapping_modules_dict = {
            method_A_vcm_id: list(set(method_A_df[method_B_vcm_id_col]))
            for method_A_vcm_id, method_A_df in method_AB_peak_overlap_df.groupby(
                method_A_vcm_id_col
            )
        }
        return method_AB_peak_overlap_df, method_AB_overlapping_modules_dict

    # def get_scores_for_methods_AB_dict(
    #         self,
    #         dataset,
    #         method_A,
    #         method_B,
    #         method_A_peaks_df,
    #         method_B_peaks_df,
    #         method_AB_peak_overlap_df,
    #         method_AB_overlapping_modules_dict,
    #         A_CMs_with_qtl=None,
    #         B_CMs_with_qtl=None,
    #         method_A_cm_id_col_in_peaks=None,
    #         method_A_peak_id_col_in_overlap=None,
    #         method_A_cm_id_col_in_overlap=None,
    #         method_B_cm_id_col_in_peaks=None,
    #         method_B_peak_id_col_in_overlap=None,
    #         method_B_cm_id_col_in_overlap=None,
    #         save_scores=False,
    #     ):
    #     if self.use_cmQTLs:
    #         if (A_CMs_with_qtl is None) and (B_CMs_with_qtl is None):
    #             sys.exit('No chromatin modules with QTLs were provided, check the input!\n Exiting...')
    #         if (A_CMs_with_qtl is not None) and (B_CMs_with_qtl is None):
    #             method_A_peaks_df = method_A_peaks_df[
    #                 method_A_peaks_df.iloc[:, method_A_cm_id_col_in_peaks].isin(A_CMs_with_qtl)
    #             ].copy()
    #         elif (A_CMs_with_qtl is None) and (B_CMs_with_qtl is not None):
    #             method_B_peaks_df = method_B_peaks_df[
    #                 method_B_peaks_df.iloc[:, method_B_cm_id_col_in_peaks].isin(B_CMs_with_qtl)
    #             ].copy()
    #         else:
    #             method_A_peaks_df = method_A_peaks_df[
    #                 method_A_peaks_df.iloc[:, method_A_cm_id_col_in_peaks].isin(A_CMs_with_qtl)
    #             ].copy()
    #             method_B_peaks_df = method_B_peaks_df[
    #                 method_B_peaks_df.iloc[:, method_B_cm_id_col_in_peaks].isin(B_CMs_with_qtl)
    #             ].copy()
    #     all_A_CMs = set(method_A_peaks_df.iloc[:, method_A_cm_id_col_in_peaks])
    #     all_B_CMs = set(method_B_peaks_df.iloc[:, method_B_cm_id_col_in_peaks])
    #     method_AB_F1_dict = {}
    #     for cm_A_id, cm_B_id_lst in tqdm.tqdm(method_AB_overlapping_modules_dict.items()):
    #         module_A_peaks = list(set(method_A_peaks_df[method_A_peaks_df.iloc[:, method_A_cm_id_col_in_peaks] == cm_A_id].iloc[:, 3]))
    #         for cm_B_id in cm_B_id_lst:
    #             module_B_peaks = list(set(method_B_peaks_df[method_B_peaks_df.iloc[:, method_B_cm_id_col_in_peaks] == cm_B_id].iloc[:, 3]))
    #             AB_overlap_df = method_AB_peak_overlap_df[
    #                 (method_AB_peak_overlap_df.iloc[:, method_A_cm_id_col_in_overlap] == cm_A_id) &
    #                 (method_AB_peak_overlap_df.iloc[:, method_B_cm_id_col_in_overlap] == cm_B_id)
    #             ]
    #             overlap_df = pd.DataFrame(
    #                 np.nan,
    #                 index=module_A_peaks,
    #                 columns=module_B_peaks
    #             )
    #             for row, col, jaccard in zip(
    #                     AB_overlap_df.iloc[:, method_A_peak_id_col_in_overlap],
    #                     AB_overlap_df.iloc[:, method_B_peak_id_col_in_overlap],
    #                     AB_overlap_df.loc[:, 'jaccard']
    #                 ):
    #                 overlap_df.loc[row, col] = jaccard
    #             overlap_df = overlap_df.fillna(0)
    #             max_in_rows = overlap_df.max(axis=1).to_numpy()
    #             max_in_cols = overlap_df.max(axis=0).to_numpy()
    #             if (self.score_type == 'HM') | (self.score_type == 'F1'):
    #                 denominator = 1 / np.mean(max_in_rows) + 1 / np.mean(max_in_cols)
    #                 if denominator != 0:
    #                     F1 = 2 / denominator
    #                 else:
    #                     F1 = 0
    #                 score = F1
    #             elif self.score_type == 'AM':
    #                 score = (np.mean(max_in_rows) + np.mean(max_in_cols)) / 2
    #             method_AB_F1_dict[cm_A_id + '_' + cm_B_id] = score
    #     A_B_overlapping_A_CMs = set(method_AB_peak_overlap_df.iloc[:, method_A_cm_id_col_in_overlap])
    #     A_B_overlapping_B_CMs = set(method_AB_peak_overlap_df.iloc[:, method_B_cm_id_col_in_overlap])
    #     method_AB_F1_dict.update(
    #         {non_overlapping_A_CMs + '_NAN': 0
    #         for non_overlapping_A_CMs in all_A_CMs - A_B_overlapping_A_CMs})
    #     method_AB_F1_dict.update({
    #         'NAN_' + non_overlapping_B_CMs: 0
    #          for non_overlapping_B_CMs in all_B_CMs - A_B_overlapping_B_CMs})
    #     if save_scores:
    #         print('Saving files with similarity scores...')
    #         overlap_output_path = os.path.join(self.output_path, self.score_type + '_scores')
    #         if not os.path.exists(overlap_output_path):
    #             os.makedirs(overlap_output_path)
    #         self.save_npy(method_A, method_B, dataset, method_AB_F1_dict, overlap_output_path)
    #     return method_AB_F1_dict

    @staticmethod
    def get_merged_peak_mapping_dict(module_peaks, new_pids):
        init_to_merged_pid_mapping = {}
        for peak_id in module_peaks:
            start, end = list(map(int, peak_id.split(":")[-2:]))
            for merged_start, merged_end in new_pids:
                if (start >= merged_start) and (end <= merged_end):
                    init_to_merged_pid_mapping[":".join([str(start), str(end)])] = (
                        ":".join([str(merged_start), str(merged_end)])
                    )
                    break
                else:
                    init_to_merged_pid_mapping[":".join([str(start), str(end)])] = (
                        ":".join([str(start), str(end)])
                    )
        return init_to_merged_pid_mapping

    def get_scores_for_methods_AB_dict(
        self,
        dataset,
        method_A,
        method_B,
        method_A_peaks_df,
        method_B_peaks_df,
        method_AB_peak_overlap_df,
        method_AB_overlapping_modules_dict,
        A_CMs_with_qtl=None,
        B_CMs_with_qtl=None,
        method_A_cm_id_col_in_peaks=None,
        method_A_peak_id_col_in_overlap=None,
        method_A_cm_id_col_in_overlap=None,
        method_B_cm_id_col_in_peaks=None,
        method_B_peak_id_col_in_overlap=None,
        method_B_cm_id_col_in_overlap=None,
        save_scores=False,
    ):
        if self.use_cmQTLs:
            if (A_CMs_with_qtl is None) and (B_CMs_with_qtl is None):
                sys.exit(
                    "No chromatin modules with QTLs were provided, check the input!\n Exiting..."
                )
            if (A_CMs_with_qtl is not None) and (B_CMs_with_qtl is None):
                method_A_peaks_df = method_A_peaks_df[
                    method_A_peaks_df.iloc[:, method_A_cm_id_col_in_peaks].isin(
                        A_CMs_with_qtl
                    )
                ].copy()
            elif (A_CMs_with_qtl is None) and (B_CMs_with_qtl is not None):
                method_B_peaks_df = method_B_peaks_df[
                    method_B_peaks_df.iloc[:, method_B_cm_id_col_in_peaks].isin(
                        B_CMs_with_qtl
                    )
                ].copy()
            else:
                method_A_peaks_df = method_A_peaks_df[
                    method_A_peaks_df.iloc[:, method_A_cm_id_col_in_peaks].isin(
                        A_CMs_with_qtl
                    )
                ].copy()
                method_B_peaks_df = method_B_peaks_df[
                    method_B_peaks_df.iloc[:, method_B_cm_id_col_in_peaks].isin(
                        B_CMs_with_qtl
                    )
                ].copy()
        all_A_CMs = set(method_A_peaks_df.iloc[:, method_A_cm_id_col_in_peaks])
        all_B_CMs = set(method_B_peaks_df.iloc[:, method_B_cm_id_col_in_peaks])
        method_AB_F1_dict = {}
        for cm_A_id, cm_B_id_lst in tqdm.tqdm(
            method_AB_overlapping_modules_dict.items()
        ):
            module_A_peaks = list(
                set(
                    method_A_peaks_df.loc[
                        method_A_peaks_df.iloc[:, method_A_cm_id_col_in_peaks]
                        == cm_A_id,
                        :,
                    ].iloc[:, 3]
                )
            )
            for cm_B_id in cm_B_id_lst:
                module_B_peaks = list(
                    set(
                        method_B_peaks_df.loc[
                            method_B_peaks_df.iloc[:, method_B_cm_id_col_in_peaks]
                            == cm_B_id,
                            :,
                        ].iloc[:, 3]
                    )
                )
                AB_overlap_df = method_AB_peak_overlap_df.loc[
                    (
                        method_AB_peak_overlap_df.iloc[:, method_A_cm_id_col_in_overlap]
                        == cm_A_id
                    )
                    & (
                        method_AB_peak_overlap_df.iloc[:, method_B_cm_id_col_in_overlap]
                        == cm_B_id
                    ),
                    :,
                ]
                merged_peak_coord_A = self.merge_intersecting_intervals(
                    [
                        list(map(int, A_peak_id.split(":")[-2:]))
                        for A_peak_id in module_A_peaks
                    ]
                )
                merged_peak_coord_B = self.merge_intersecting_intervals(
                    [
                        list(map(int, B_peak_id.split(":")[-2:]))
                        for B_peak_id in module_B_peaks
                    ]
                )
                new_pids_A = [
                    ":".join([str(A_start), str(A_end)])
                    for A_start, A_end in merged_peak_coord_A
                ]
                new_pids_B = [
                    ":".join([str(B_start), str(B_end)])
                    for B_start, B_end in merged_peak_coord_B
                ]
                init_to_merged_pid_mapping_A = self.get_merged_peak_mapping_dict(
                    module_A_peaks, merged_peak_coord_A
                )
                init_to_merged_pid_mapping_B = self.get_merged_peak_mapping_dict(
                    module_B_peaks, merged_peak_coord_B
                )

                overlap_df = pd.DataFrame(0, index=new_pids_A, columns=new_pids_B)
                for row, col, jaccard in zip(
                    AB_overlap_df.iloc[:, method_A_peak_id_col_in_overlap],
                    AB_overlap_df.iloc[:, method_B_peak_id_col_in_overlap],
                    AB_overlap_df.loc[:, "jaccard"],
                ):
                    row_peak_str = ":".join(row.split(":")[-2:])
                    col_peak_str = ":".join(col.split(":")[-2:])
                    row = init_to_merged_pid_mapping_A[row_peak_str]
                    col = init_to_merged_pid_mapping_B[col_peak_str]
                    overlap_df.loc[row, col] += jaccard
                overlap_df[overlap_df > 1] = 1
                max_in_rows = overlap_df.max(axis=1).to_numpy()
                max_in_cols = overlap_df.max(axis=0).to_numpy()
                if (self.score_type == "HM") | (self.score_type == "F1"):
                    denominator = 1 / np.mean(max_in_rows) + 1 / np.mean(max_in_cols)
                    if denominator != 0:
                        F1 = 2 / denominator
                    else:
                        F1 = 0
                    score = F1
                elif self.score_type == "AM":
                    score = (np.mean(max_in_rows) + np.mean(max_in_cols)) / 2
                method_AB_F1_dict[cm_A_id + "_" + cm_B_id] = score
        A_B_overlapping_A_CMs = set(
            method_AB_peak_overlap_df.iloc[:, method_A_cm_id_col_in_overlap]
        )
        A_B_overlapping_B_CMs = set(
            method_AB_peak_overlap_df.iloc[:, method_B_cm_id_col_in_overlap]
        )
        method_AB_F1_dict.update(
            {
                non_overlapping_A_CMs + "_NAN": 0
                for non_overlapping_A_CMs in all_A_CMs - A_B_overlapping_A_CMs
            }
        )
        method_AB_F1_dict.update(
            {
                "NAN_" + non_overlapping_B_CMs: 0
                for non_overlapping_B_CMs in all_B_CMs - A_B_overlapping_B_CMs
            }
        )
        if save_scores:
            print("Saving files with similarity scores...")
            overlap_output_path = os.path.join(
                self.output_path, self.score_type + "_scores"
            )
            if not os.path.exists(overlap_output_path):
                os.makedirs(overlap_output_path)
            self.save_npy(
                method_A, method_B, dataset, method_AB_F1_dict, overlap_output_path
            )
        return method_AB_F1_dict

    def subset_tracks_df(
        self,
        method,
        dataset,
        CMs_to_subset,
        similarity_threshold,
        path_to_cm_tracks_method,
        min_cm_size=0,
        save_filtered_tracks=False,
        direction=None,
        zero_cms_to_exclude=None,
    ):
        tracks_df = pd.read_csv(path_to_cm_tracks_method, sep="\t", header=None)
        if min_cm_size != 0:
            tracks_df = tracks_df[tracks_df.iloc[:, 9] >= min_cm_size]
        tracks_df.iloc[:, 3] = tracks_df.iloc[:, 3].str.replace("_", "~")
        if (direction == "less") or (direction == "lower"):
            method_specific_CMs_df_le = tracks_df.loc[
                ~tracks_df.iloc[:, 3].isin(CMs_to_subset)
            ]
            method_specific_CMs_df_le = method_specific_CMs_df_le.loc[
                ~tracks_df.iloc[:, 3].isin(zero_cms_to_exclude), :
            ].copy()
            method_specific_CMs_df_gr = None
        elif (direction == "larger") or (direction == "greater"):
            method_specific_CMs_df_le = None
            method_specific_CMs_df_gr = tracks_df[
                tracks_df.iloc[:, 3].isin(CMs_to_subset)
            ]
        elif direction == "both":
            method_specific_CMs_df_le = tracks_df[
                ~tracks_df.iloc[:, 3].isin(CMs_to_subset)
            ]
            method_specific_CMs_df_gr = tracks_df[
                tracks_df.iloc[:, 3].isin(CMs_to_subset)
            ]
        else:
            direction = None
            method_specific_CMs_df_le = None
            method_specific_CMs_df_gr = tracks_df[
                tracks_df.iloc[:, 3].isin(CMs_to_subset)
            ]

        if save_filtered_tracks:
            if self.use_cmQTLs:
                folder_name = "modules_with_QTL"
            else:
                folder_name = "all_modules"

            tracks_output_path = os.path.join(
                self.output_path, "subset_tracks_and_content", folder_name
            )
            if not os.path.exists(tracks_output_path):
                os.makedirs(tracks_output_path)
            if direction is not None:
                if method_specific_CMs_df_le is not None:
                    self.save_tracks(
                        method_specific_CMs_df_le,
                        dataset,
                        method,
                        similarity_threshold,
                        tracks_output_path,
                        direction="less",
                    )
                if method_specific_CMs_df_gr is not None:
                    self.save_tracks(
                        method_specific_CMs_df_gr,
                        dataset,
                        method,
                        similarity_threshold,
                        tracks_output_path,
                        direction="greater",
                    )
            else:
                self.save_tracks(
                    method_specific_CMs_df_gr,
                    dataset,
                    method,
                    similarity_threshold,
                    tracks_output_path,
                    direction=None,
                )
        return method_specific_CMs_df_le, method_specific_CMs_df_gr

    def get_CM_similarity_scores_df(
        self,
        dataset,
        similarity_threshold,
        method_A,
        method_B,
        A_B_F1_scores,
        save_df=False,
        unique=False,
    ):
        CMs_and_similarity_scores = []
        for A_B_pair, A_B_F1 in A_B_F1_scores.items():
            cm_id_A, cm_id_B = A_B_pair.split("_")[0], "_".join(A_B_pair.split("_")[1:])
            CMs_and_similarity_scores.append(
                [method_A + "_" + cm_id_A, method_B + "_" + cm_id_B, A_B_F1]
            )
        AB_CMs_similarity_score_df = pd.DataFrame(
            CMs_and_similarity_scores, columns=["module_1", "module_2", "score"]
        )
        # Select the most universal chromatin modules
        if unique:
            return (
                AB_CMs_similarity_score_df[
                    AB_CMs_similarity_score_df.loc[:, "score"] > similarity_threshold
                ],
                None,
            )
        if (similarity_threshold > 0) and (similarity_threshold < 1):
            suffix = "greater"
            filt_AB_CMs_similarity_score_df = AB_CMs_similarity_score_df[
                AB_CMs_similarity_score_df.loc[:, "score"] > similarity_threshold
            ]
            # elif suffix == 'less':
            #     filt_AB_CMs_similarity_score_df = AB_CMs_similarity_score_df[
            #         (AB_CMs_similarity_score_df.loc[:, 'score'] < similarity_threshold) &
            #         (AB_CMs_similarity_score_df.loc[:, 'score'] != 0)
            #     ]
            zero_similarity_score_df = AB_CMs_similarity_score_df[
                AB_CMs_similarity_score_df.loc[:, "score"] == 0
            ]
        elif (similarity_threshold == 0) or (similarity_threshold == 1):
            suffix = "equal"
            filt_AB_CMs_similarity_score_df = AB_CMs_similarity_score_df[
                AB_CMs_similarity_score_df.loc[:, "score"] == similarity_threshold
            ]
            zero_similarity_score_df = None
        elif similarity_threshold is None:
            suffix = "all"
            similarity_threshold = ""
            filt_AB_CMs_similarity_score_df = AB_CMs_similarity_score_df
            zero_similarity_score_df = None
        if save_df:
            output_file_name = os.path.join(
                self.output_path,
                "_".join(
                    [
                        dataset,
                        method_A,
                        method_B,
                        self.score_type,
                        "scores_for_CM_pairs",
                        suffix,
                        str(similarity_threshold) + ".bed",
                    ]
                ),
            )
            filt_AB_CMs_similarity_score_df.to_csv(
                output_file_name, sep="\t", header=False, index=False
            )
        return filt_AB_CMs_similarity_score_df, zero_similarity_score_df

    @staticmethod
    def get_max_cliques(filt_CMs_similarity_score_df, clique_size):
        G = nx.from_pandas_edgelist(
            filt_CMs_similarity_score_df,
            source="module_1",
            target="module_2",
            edge_attr="score",
            create_using=nx.Graph,
        )
        CMs_to_subset = []
        for clique_nodes in nx.find_cliques(G):
            g = G.subgraph(clique_nodes)
            if g.number_of_nodes() == clique_size:
                CMs_to_subset.extend(clique_nodes)
        return CMs_to_subset

    def get_CM_similarity_score_df_for_all_methods(
        self, dataset, similarity_threshold, unique
    ):
        all_method_pairs = list(itertools.permutations(self.methods_list, 2))
        all_df_list = []
        all_zero_cms_list = []
        for method_A, method_B in all_method_pairs:
            A_B_F1_scores = self.load_npy(
                method_A,
                method_B,
                dataset,
                os.path.join(self.output_path, self.score_type + "_scores"),
            )
            (
                filt_AB_CMs_similarity_score_df,
                zero_cms_df,
            ) = self.get_CM_similarity_scores_df(
                dataset,
                similarity_threshold,
                method_A,
                method_B,
                A_B_F1_scores,
                unique=unique,
                save_df=False,
            )
            all_df_list.append(filt_AB_CMs_similarity_score_df)
            all_zero_cms_list.append(zero_cms_df)
        if unique or (similarity_threshold == 0) or (similarity_threshold == 1):
            return pd.concat(all_df_list), None
        else:
            return pd.concat(all_df_list), pd.concat(all_zero_cms_list)

    def get_unique_CMs_df(
        self,
        method,
        dataset,
        cm_tracks_path,
        similarity_threshold,
        min_cm_size=0,
        save_filtered_tracks=False,
    ):
        (
            all_ds_CMs_similarity_scores_greater_zero,
            _,
        ) = self.get_CM_similarity_score_df_for_all_methods(
            dataset, similarity_threshold=0, unique=True
        )
        (
            filt_all_ds_CMs_similarity_score_df,
            _,
        ) = self.get_CM_similarity_score_df_for_all_methods(
            dataset, similarity_threshold=0, unique=False
        )
        CMs_to_exclude = list(
            set(list(all_ds_CMs_similarity_scores_greater_zero.iloc[:, 0])).union(
                set(list(all_ds_CMs_similarity_scores_greater_zero.iloc[:, 1]))
            )
        )
        CMs_to_subset = list(
            set(list(filt_all_ds_CMs_similarity_score_df.iloc[:, 0])).union(
                set(list(filt_all_ds_CMs_similarity_score_df.iloc[:, 1]))
            )
        )
        method_specific_CMs_to_subset = [
            cm_to_subset.split("_")[1]
            for cm_to_subset in CMs_to_subset
            if (method in cm_to_subset) and not (cm_to_subset in CMs_to_exclude)
        ]
        _, method_specific_CMs_df_gr = self.subset_tracks_df(
            method,
            dataset,
            method_specific_CMs_to_subset,
            similarity_threshold,
            cm_tracks_path,
            save_filtered_tracks=save_filtered_tracks,
            min_cm_size=min_cm_size,
            direction=None,
        )
        return _, method_specific_CMs_df_gr

    def get_universal_CMs_df(
        self,
        method,
        dataset,
        cm_tracks_path,
        similarity_threshold,
        min_cm_size=0,
        save_filtered_tracks=False,
        direction=None,
    ):
        (
            filt_all_ds_CMs_similarity_score_df,
            zero_cms,
        ) = self.get_CM_similarity_score_df_for_all_methods(
            dataset, similarity_threshold, unique=False
        )

        zero_cms_to_exclude = None
        if similarity_threshold > 0:
            CMs_to_subset = self.get_max_cliques(
                filt_all_ds_CMs_similarity_score_df, clique_size=self.n_datasets
            )
            if zero_cms is None:
                zero_cms_to_exclude = []
            else:
                zero_cms_to_exclude = [
                    zero_cm_id
                    for zero_pair in list(zip(zero_cms.iloc[:, 0], zero_cms.iloc[:, 1]))
                    for zero_cm_id in zero_pair
                ]
        else:
            sys.exit("Similarity threshold in not greater than zero!")
        method_specific_CMs_to_subset = [
            cm_to_subset.split("_")[1]
            for cm_to_subset in CMs_to_subset
            if method in cm_to_subset
        ]
        method_specific_CMs_df_le, method_specific_CMs_df_gr = self.subset_tracks_df(
            method,
            dataset,
            method_specific_CMs_to_subset,
            similarity_threshold,
            cm_tracks_path,
            save_filtered_tracks=save_filtered_tracks,
            direction=direction,
            min_cm_size=min_cm_size,
            zero_cms_to_exclude=zero_cms_to_exclude,
        )
        return method_specific_CMs_df_le, method_specific_CMs_df_gr

    # @staticmethod
    # def get_tracks_path(cm_method, cell_type, path_to_cms, pp_threshold=0.8):
    #     if cm_method == "clomics":
    #         method_name = "Clomics"
    #         tracks_path = os.path.join(path_to_cms, "vcm.tracks.bed")
    #     elif cm_method == "vcmtools":
    #         vcm_pv_threshold = 0.001
    #         method_name = "VCMtools"
    #         tracks_path = os.path.join(
    #             path_to_cms,
    #             "_".join(
    #                 [
    #                     cell_type,
    #                     "NOT_TOTEM_VCMs_corrected_pvalue",
    #                     str(vcm_pv_threshold) + ".tracks.bed",
    #                 ]
    #             ),
    #         )
    #     elif cm_method == "PHM":
    #         method_name = "PHM"
    #         tracks_path = os.path.join(
    #             path_to_cms,
    #             "_".join(
    #                 [
    #                     cell_type,
    #                     str(pp_threshold),
    #                     "merged_phm_all_chr_NOT_TOTEM.tracks.bed",
    #                 ]
    #             ),
    #         )
    #     return method_name, tracks_path

    # @staticmethod
    # def get_genes_overlapping_CMs(
    #         CM_gene_dict_path,
    #         peak_to_cm_mapping,
    #         subset_CMs=False,
    #         hgnc=False,
    #         CM_list=None,
    #         gene_transcript_mapping_df=None
    #     ):
    #     gene_body_and_promoter_dict = np.load(
    #         CM_gene_dict_path,
    #         allow_pickle=True
    #     ).item()
    #     if subset_CMs:
    #         genes_overlapping_modules = []
    #         for pid, gene_list in gene_body_and_promoter_dict.items():
    #             cm_id = peak_to_cm_mapping.get(pid)
    #             if cm_id in CM_list:
    #                 genes_overlapping_modules.extend(gene_list)
    #     else:
    #         genes_overlapping_modules = [
    #             gene_id
    #             for gene_lst in gene_body_and_promoter_dict.values()
    #             for gene_id in gene_lst
    #         ]
    #     if hgnc:
    #         hgnc_genes_overlapping_modules = gene_transcript_mapping_df[
    #             gene_transcript_mapping_df['ensembl_gene_id'].isin(genes_overlapping_modules)
    #             ]['hgnc_symbol'].copy()
    #         return list(set([
    #             hgnc_id
    #             for hgnc_id in hgnc_genes_overlapping_modules
    #             if str(hgnc_id) != 'nan'
    #             ]))
    #     else:
    #         return list(set(genes_overlapping_modules))

    # def get_enriched_terms(
    #         self,
    #         method,
    #         dataset,
    #         similarity_threshold,
    #         genes_overlapping_modules,
    #         enrichr_libraries,
    #         process_libraries_together=False,
    #         pv_filtered_output=False,
    #         enrichment_pv_threshold=0.05,
    #         save_enrichr_df=False,
    #         how=None,
    #         direction=None
    #     ):
    #     if not genes_overlapping_modules:
    #         print('There are no genes overlapping chromatin modules')
    #         return None
    #     if process_libraries_together:
    #         print('Running Enrichr on all libraries...')
    #         enr = gp.enrichr(
    #                 gene_list=genes_overlapping_modules,
    #                 gene_sets=enrichr_libraries,
    #                 organism='Human',
    #                 description='_'.join([dataset, method, 'enrichr']),
    #                 outdir=os.path.join(
    #                     self.output_path,
    #                     'temp',
    #                     '_'.join([dataset, method, 'enrichr']))
    #             )
    #         if enr.results.shape[0] == 0:
    #                 print('There are no enriched terms for the gene set. Stopping...')
    #                 return None
    #         if pv_filtered_output:
    #             enrichr_results_df = enr.results
    #             enrichr_results_df_filt = enrichr_results_df[
    #                 enrichr_results_df['Adjusted P-value'] <= enrichment_pv_threshold]
    #             enrichr_results_df_filt.loc[:, '-log10(Adjusted pv)'] = -np.log10(enrichr_results_df_filt['Adjusted P-value'])
    #             enrichr_results_df_filt = enrichr_results_df_filt[['Term', '-log10(Adjusted pv)', '']]
    #         else:
    #             enrichr_results_df_filt = enr.results
    #         if (save_enrichr_df) and (enrichr_results_df_filt.shape[0] != 0):
    #             self.save_enrichr_df(
    #                 method,
    #                 dataset,
    #                 similarity_threshold,
    #                 enrichr_results_df_filt,
    #                 'all_input_libraries_together',
    #                 pv_filtered_output,
    #                 enrichment_pv_threshold,
    #                 how=how,
    #                 direction=direction
    #                 )
    #         return enrichr_results_df_filt
    #     else:
    #         print('Running Enrichr for each library separately...')
    #         enrichr_results_per_library = {}
    #         for enrichr_library in enrichr_libraries:
    #             enr = gp.enrichr(
    #                     gene_list=genes_overlapping_modules,
    #                     gene_sets=enrichr_library,
    #                     organism='Human',
    #                     description='_'.join([dataset, method, 'enrichr']),
    #                     outdir=os.path.join(
    #                         self.output_path,
    #                         'temp',
    #                         '_'.join([dataset, method, 'enrichr']))
    #                 )
    #             if enr.results.shape[0] == 0:
    #                 print('There are no enriched terms for the gene set. Stopping...')
    #                 return None
    #             if pv_filtered_output:
    #                 enrichr_results_df = enr.results
    #                 enrichr_results_df_filt = enrichr_results_df[enrichr_results_df['Adjusted P-value'] <= enrichment_pv_threshold]
    #             else:
    #                 enrichr_results_df_filt = enr.results
    #             enrichr_results_df_filt.loc[:, '-log10(Adjusted pv)'] = -np.log10(enrichr_results_df_filt['Adjusted P-value'])
    #             enrichr_results_df_filt.loc[:, 'n_hits'] = enrichr_results_df_filt['Genes'].str.split(';').str.len()
    #             if (save_enrichr_df) and (enrichr_results_df_filt.shape[0] != 0):
    #                 self.save_enrichr_df(
    #                     method,
    #                     dataset,
    #                     similarity_threshold,
    #                     enrichr_results_df_filt,
    #                     enrichr_library,
    #                     pv_filtered_output,
    #                     enrichment_pv_threshold,
    #                     how=how,
    #                     direction=direction
    #                     )
    #                 gene_counts_df = pd.DataFrame.from_dict(collections.Counter([
    #                     gene
    #                     for gene_str in list(enrichr_results_df_filt['Genes'])
    #                     for gene in gene_str.split(';')
    #                 ]),
    #                 orient='index').reset_index()
    #                 gene_counts_df.columns = ['hgnc_symbol', 'occurence']
    #                 self.save_enriched_genes(
    #                     method,
    #                     dataset,
    #                     similarity_threshold,
    #                     gene_counts_df,
    #                     enrichr_library,
    #                     pv_filtered_output,
    #                     enrichment_pv_threshold,
    #                     direction=direction
    #                 )
    #             enrichr_results_per_library[enrichr_library] = enrichr_results_df_filt
    #         return enrichr_results_per_library
