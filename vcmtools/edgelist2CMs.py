import numpy as np
import os
import pandas as pd
import statsmodels.stats.multitest as multi
import tqdm


class edgelist2CMs:
    def __init__(self, output_path, remove_totem, save_files, save_meta):
        self.output_path = output_path
        self.remove_totem = remove_totem
        self.save_files = save_files
        self.save_meta = save_meta

    @staticmethod
    def connected_tuples(pairs):
        # for every element, we keep a reference to the list it belongs to
        lists_by_element = {}

        def make_new_list_for(x, y):
            lists_by_element[x] = lists_by_element[y] = [x, y]

        def add_element_to_list(lst, el):
            lst.append(el)
            lists_by_element[el] = lst

        def merge_lists(lst1, lst2):
            merged_list = lst1 + lst2
            for el in merged_list:
                lists_by_element[el] = merged_list

        for x, y in pairs:
            xList = lists_by_element.get(x)
            yList = lists_by_element.get(y)

            if not xList and not yList:
                make_new_list_for(x, y)

            if xList and not yList:
                add_element_to_list(xList, y)

            if yList and not xList:
                add_element_to_list(yList, x)

            if xList and yList and xList != yList:
                merge_lists(xList, yList)
        return set(tuple(l) for l in lists_by_element.values())

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

    def save_dict(self, file_name, dataset, dictionary):
        if not os.path.exists(os.path.join(self.output_path, dataset)):
            os.makedirs(os.path.join(self.output_path, dataset))
        np.save(os.path.join(self.output_path, dataset, file_name), dictionary)

    @staticmethod
    def load_corr_edgelist(
        path_to_edgelist,
    ):
        all_corr = pd.read_csv(path_to_edgelist, sep="\t", header=None)
        all_corr.columns = ["peak1", "peak2", "corr", "p_value"]
        return all_corr

    @staticmethod
    def correct_pv(init_pv_df, pv_threshold):
        all_pv_df = init_pv_df.copy()
        _, pvals_corrected, _, _ = multi.multipletests(
            all_pv_df.loc[:, "p_value"],
            alpha=pv_threshold,
            method="fdr_bh",  # BH correction
        )
        all_pv_df.loc[:, "pvals_corrected"] = pvals_corrected
        return all_pv_df[all_pv_df.loc[:, "pvals_corrected"] <= pv_threshold]

    def get_CMs(self, dataset, pv_filtered_corr_edgelist, pvalue_threshold):
        connected_components = self.connected_tuples(
            list(zip(pv_filtered_corr_edgelist.peak1, pv_filtered_corr_edgelist.peak2))
        )
        bed_df = pd.DataFrame(
            columns=[
                "chr",
                "start",
                "end",
                "CM_id",
                "number",
                "strain",
                "start_duplicate",
                "end_duplicate",
                "numbers",
                "CM_size",
                "peak_length",
                "peak_starts",
                "totem_CM",
            ],
            index=np.arange(len(connected_components)),
        )
        cm_content = pd.DataFrame(
            columns=["CM_id", "CM_size", "peaks"],
            index=np.arange(len(connected_components)),
        )
        for index, cm in tqdm.tqdm(enumerate(connected_components)):
            start_end = sorted(
                [(int(peak.split(":")[-2]), int(peak.split(":")[-1])) for peak in cm],
                key=lambda x: x[0],
            )

            start = np.array([pair[0] for pair in start_end])
            peak_starts = list(start - min(start_end, key=lambda x: x[0])[0])

            cm_content["CM_id"].loc[index] = "cm" + str(index)
            cm_content["CM_size"].loc[index] = int(len(cm))
            cm_content["peaks"].loc[index] = ",".join(cm)

            if "chr" in str(cm[0].split(":")[1]):
                bed_df["chr"].loc[index] = str(cm[0].split(":")[1])
            else:
                bed_df["chr"].loc[index] = "chr" + str(cm[0].split(":")[1])
            bed_df["start"].loc[index] = min(start_end, key=lambda x: x[0])[0]
            bed_df["end"].loc[index] = max(start_end, key=lambda x: x[1])[1]
            bed_df["CM_id"].loc[index] = "cm" + str(index)
            bed_df["number"].loc[index] = int(1000)
            bed_df["strain"].loc[index] = "+"
            bed_df["start_duplicate"].loc[index] = min(start_end, key=lambda x: x[0])[0]
            bed_df["end_duplicate"].loc[index] = max(start_end, key=lambda x: x[1])[1]
            bed_df["numbers"].loc[index] = "0,0,0"
            bed_df["CM_size"].loc[index] = int(len(cm))
            bed_df["peak_length"].loc[index] = ",".join(
                map(str, [j - i for i, j in start_end])
            )
            bed_df["peak_starts"].loc[index] = ",".join(map(str, peak_starts))
            bed_df["totem_CM"].loc[index] = (
                1 if len(self.merge_intersecting_intervals(start_end)) == 1 else 0
            )

        if bed_df.shape[0] == 0:
            return pd.DataFrame(), pd.DataFrame(), {}

        if not "chr" in cm_content.iloc[0, 2]:
            for i in list(range(1, 23)) + ["X", "Y"]:
                cm_content["peaks"] = cm_content["peaks"].str.replace(
                    ":" + str(i) + ":", ":chr" + str(i) + ":"
                )
        sorted_bed_df = bed_df.sort_values(["chr", "start"]).copy()
        cm_content = cm_content.set_index("CM_id")
        sorted_cm_content = cm_content.loc[sorted_bed_df["CM_id"].copy()].copy()
        sorted_bed_df["CM_id"] = [
            "cm" + str(cm_int_idx) for cm_int_idx in np.arange(sorted_bed_df.shape[0])
        ]
        sorted_cm_content["new_CM_id"] = [
            "cm" + str(cm_int_idx)
            for cm_int_idx in np.arange(sorted_cm_content.shape[0])
        ]
        sorted_cm_content = sorted_cm_content[["new_CM_id", "CM_size", "peaks"]]
        sorted_cm_content.columns = ["CM_id", "CM_size", "peaks"]
        if self.remove_totem:
            final_cm_tracks = sorted_bed_df[sorted_bed_df.loc[:, "totem_CM"] == 0]
            final_cm_content = sorted_cm_content[
                sorted_cm_content.loc[:, "CM_id"].isin(final_cm_tracks.loc[:, "CM_id"])
            ]
        else:
            final_cm_tracks = sorted_bed_df
            final_cm_content = sorted_cm_content
        meta_dict = dict()
        if self.save_meta:
            if len(connected_components) != 0:
                meta_dict["nCMs"] = len(connected_components)
                meta_dict["mean_abs_corr"] = np.mean(
                    abs(pv_filtered_corr_edgelist["corr"])
                )
                meta_dict["CM_lengths"] = list(
                    final_cm_tracks.loc[:, "end"] - final_cm_tracks.loc[:, "start"]
                )
                meta_dict["CM_sizes"] = list(final_cm_tracks.loc[:, "CM_size"])
            else:
                meta_dict["nCMs"] = 0
                meta_dict["mean_abs_corr"] = np.nan
                meta_dict["CM_lengths"] = []
                meta_dict["CM_sizes"] = []

        if self.save_files:
            if not os.path.exists(os.path.join(self.output_path, dataset)):
                os.makedirs(os.path.join(self.output_path, dataset))
            if self.remove_totem:
                suffix = "NOT_TOTEM"
                sorted_bed_df.to_csv(
                    os.path.join(
                        self.output_path,
                        dataset,
                        "_".join(
                            [
                                dataset,
                                "all_VCMtools_CMs_corrected_pvalue",
                                str(pvalue_threshold) + ".tracks.bed",
                            ]
                        ),
                    ),
                    sep="\t",
                    header=False,
                    index=False,
                )
                sorted_cm_content.to_csv(
                    os.path.join(
                        self.output_path,
                        dataset,
                        "_".join(
                            [
                                dataset,
                                "all_VCMtools_CMs_corrected_pvalue",
                                str(pvalue_threshold) + ".content.txt",
                            ]
                        ),
                    ),
                    sep="\t",
                    header=False,
                    index=False,
                )
            else:
                suffix = "all"
            final_cm_tracks.to_csv(
                os.path.join(
                    self.output_path,
                    dataset,
                    "_".join(
                        [
                            dataset,
                            suffix,
                            "VCMtools_CMs_corrected_pvalue",
                            str(pvalue_threshold) + ".tracks.bed",
                        ]
                    ),
                ),
                sep="\t",
                header=False,
                index=False,
            )
            final_cm_content.to_csv(
                os.path.join(
                    self.output_path,
                    dataset,
                    "_".join(
                        [
                            dataset,
                            suffix,
                            "VCMtools_CMs_corrected_pvalue",
                            str(pvalue_threshold) + ".content.txt",
                        ]
                    ),
                ),
                sep="\t",
                header=False,
                index=False,
            )
        return final_cm_tracks, final_cm_content, meta_dict
