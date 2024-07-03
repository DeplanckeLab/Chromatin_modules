import numpy as np
import pandas as pd
import rpy2
import sys
from collections import Counter
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()


class Simulation:
    def __init__(
        self,
        method,
        n_sim_per_cm,
        chromosome,
        cm_df,
        cm_peak_coord_dict,
        ncm_peak_coord_dict,
        ncm_k27ac_peak_coord_dict,
        ncm_k4me1_peak_coord_dict,
        chr_count_mtx,
        single_lead_peaks=None,
    ):
        self.method = method
        self.n_sim_per_cm = n_sim_per_cm
        self.chromosome = chromosome
        self.cm_df = cm_df.reset_index()
        self.cm_peak_coord_dict = cm_peak_coord_dict
        self.ncm_peak_coord_dict = ncm_peak_coord_dict
        self.ncm_k27ac_peak_coord_dict = ncm_k27ac_peak_coord_dict
        self.ncm_k4me1_peak_coord_dict = ncm_k4me1_peak_coord_dict
        self.chr_count_mtx = chr_count_mtx
        if self.method in ["PHM", "phm"]:
            self.single_lead_peaks = single_lead_peaks

    # Function to find maximal disjoint set
    @staticmethod
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

    def get_simulated_cm_peak_ids(
        self, all_peaks, set_of_used_peaks, window, n_peaks_in_cm, mark_peak_dict, mark
    ):
        if self.method in ["PHM", "phm"]:
            chromosome_str = str(self.chromosome)
            if not "chr" in chromosome_str:
                chromosome_str = ":chr" + chromosome_str + ":"
            else:
                chromosome_str = ":" + chromosome_str + ":"
            peak_id = np.random.choice(
                [
                    single_peak
                    for single_peak in list(self.single_lead_peaks)
                    if pd.Series([single_peak]).str.contains(mark)[0]
                    & (chromosome_str in single_peak)
                ],
                1,
                replace=False,
            )[0]
        else:
            peak_id = np.random.choice(list(all_peaks), 1, replace=False)[0]
        peak_coord = mark_peak_dict[peak_id]
        proximal_peaks = [
            peak_2
            for peak_2 in all_peaks
            if (
                (peak_coord[0] + peak_coord[1]) / 2 - window
                <= mark_peak_dict[peak_2][0]
            )
            & (
                mark_peak_dict[peak_2][1]
                <= (peak_coord[0] + peak_coord[1]) / 2 + window
            )
            & (peak_2 not in set_of_used_peaks)
        ]
        if (not proximal_peaks) or (len(proximal_peaks) < n_peaks_in_cm):
            return list()
        if len(proximal_peaks) == n_peaks_in_cm:
            simulated_cm_peak_ids = proximal_peaks
        else:
            proximal_peak_coord = [(pid, mark_peak_dict[pid]) for pid in proximal_peaks]
            proximal_peak_coord.sort(key=lambda x: x[1])
            (first_pid, _), (last_pid, _) = (
                proximal_peak_coord[0],
                proximal_peak_coord[-1],
            )
            remaining_peaks = list(
                set(proximal_peaks) - set([first_pid]) - set([last_pid])
            )
            sampled_proximal_peaks = list(
                set(np.random.choice(remaining_peaks, n_peaks_in_cm - 2, replace=False))
            )
            simulated_cm_peak_ids = list(
                set([first_pid, last_pid] + sampled_proximal_peaks)
            )
        return simulated_cm_peak_ids

    @staticmethod
    def get_cm_gc_content(peak_id_list):
        bs_genome = importr("BSgenome.Hsapiens.UCSC.hg19")
        bio_strings = importr("Biostrings")
        genomic_ranges = importr("GenomicRanges")
        peak_id_df = pd.DataFrame(
            [
                [p_id.split(":")[0], int(p_id.split(":")[1]), int(p_id.split(":")[2])]
                for p_id in peak_id_list
            ]
        )
        peak_id_df.columns = ["seqnames", "start", "end"]
        peak_id_df.index = peak_id_list
        peak_id_df = peak_id_df.drop_duplicates()
        rpy2.robjects.globalenv["peaks_df"] = peak_id_df
        granges_for_peaks = genomic_ranges.makeGRangesFromDataFrame(peak_id_df)
        sequences_for_peaks = bio_strings.getSeq(bs_genome.Hsapiens, granges_for_peaks)
        rpy2.robjects.globalenv["sequences_for_peaks"] = sequences_for_peaks
        gc_content = rpy2.robjects.r(
            'as.numeric(Biostrings::letterFrequency(x = sequences_for_peaks, letters = "GC", as.prob = TRUE))'
        )
        return np.mean(gc_content)

    def get_simulated_cm_peak_ids_both_marks(
        self,
        k27ac_peaks,
        k4me1_peaks,
        number_of_k27ac_peaks,
        set_of_used_peaks,
        window,
        n_peaks_in_cm,
    ):
        if self.method == "PHM":
            k27ac_peak = np.random.choice(
                list(set(k27ac_peaks).intersection(set(self.single_lead_peaks))),
                1,
                replace=False,
            )[0]
        else:
            k27ac_peak = np.random.choice(list(k27ac_peaks), 1, replace=False)[0]

        k27ac_peak_coord = self.ncm_k27ac_peak_coord_dict[k27ac_peak]
        proximal_k27ac_peaks = [
            k27ac_peak_2
            for k27ac_peak_2 in k27ac_peaks
            if (
                (k27ac_peak_coord[0] + k27ac_peak_coord[1]) / 2 - window
                <= self.ncm_k27ac_peak_coord_dict[k27ac_peak_2][0]
            )
            & (
                self.ncm_k27ac_peak_coord_dict[k27ac_peak_2][1]
                <= (k27ac_peak_coord[0] + k27ac_peak_coord[1]) / 2 + window
            )
            & (k27ac_peak_2 not in set_of_used_peaks)
        ]
        proximal_k4me1_peaks = [
            k4me1_peak_2
            for k4me1_peak_2 in k4me1_peaks
            if (
                (k27ac_peak_coord[0] + k27ac_peak_coord[1]) / 2 - window
                <= self.ncm_k4me1_peak_coord_dict[k4me1_peak_2][0]
            )
            & (
                self.ncm_k4me1_peak_coord_dict[k4me1_peak_2][1]
                <= (k27ac_peak_coord[0] + k27ac_peak_coord[1]) / 2 + window
            )
            & (k4me1_peak_2 not in set_of_used_peaks)
        ]
        if len(proximal_k27ac_peaks) < number_of_k27ac_peaks:
            return list()
        if len(proximal_k4me1_peaks) < n_peaks_in_cm - number_of_k27ac_peaks:
            return list()
        sampled_k27ac_proximal_peaks = list(
            set(
                np.random.choice(
                    proximal_k27ac_peaks, number_of_k27ac_peaks, replace=False
                )
            )
        )
        sampled_k4me1_proximal_peaks = list(
            set(
                np.random.choice(
                    proximal_k4me1_peaks,
                    n_peaks_in_cm - number_of_k27ac_peaks,
                    replace=False,
                )
            )
        )
        simulated_cm_peak_ids = list(
            set(sampled_k27ac_proximal_peaks + sampled_k4me1_proximal_peaks)
        )
        return simulated_cm_peak_ids

    def run_one_simulation(self):
        sampled_cm_parameters = list(
            zip(
                self.cm_df["cm_id"],
                self.cm_df["cm_size"],
                self.cm_df["cm_length"],
                self.cm_df["H3K27ac_fraction"],
            )
        )
        all_k27ac_peaks = set(self.ncm_k27ac_peak_coord_dict.keys())
        all_k4me1_peaks = set(self.ncm_k4me1_peak_coord_dict.keys())

        used_peaks = set()
        set_of_all_peaks = all_k27ac_peaks.union(all_k4me1_peaks)
        ref_feature_vec = {}
        simulated_cms_pid_dict = {}
        simulated_cms_dict = {}
        sim_feature_vec_dict = {}
        for list_of_params in sampled_cm_parameters:
            cm_id, cm_size, cm_length, ref_k27ac_fraction = list_of_params
            eps = cm_length + 2500  # abs(np.random.normal()) * 500

            ref_cm_peaks = sorted(self.cm_peak_coord_dict[cm_id])
            ref_cm_peak_bp = np.sum([p[1] - p[0] for p in ref_cm_peaks])
            ref_pid = self.get_pid(ref_cm_peaks)
            ref_mean_mean, _, ref_mean_std, _ = self.get_peak_stats(ref_pid)
            ref_feature_dict = {
                "cm_size": cm_size,
                "cm_length": cm_length,
                "peak_mean": ref_mean_mean,
                "peak_std": ref_mean_std,
                "peak_bp_length": ref_cm_peak_bp,
                "k27ac_fraction": ref_k27ac_fraction,
            }

            n_k27ac_peaks = round(ref_k27ac_fraction * cm_size)

            ref_feature_vec[cm_id] = ref_feature_dict
            simulated_cms_pids = {}
            simulated_cms = {}
            sim_feature_vec = {}
            for i in np.arange(self.n_sim_per_cm):
                if n_k27ac_peaks == 0:
                    simulated_cm_peak_ids = self.get_simulated_cm_peak_ids(
                        all_k4me1_peaks,
                        used_peaks,
                        eps,
                        cm_size,
                        self.ncm_k4me1_peak_coord_dict,
                        mark="H3K4me1|ATAC",
                    )
                elif n_k27ac_peaks == cm_size:
                    simulated_cm_peak_ids = self.get_simulated_cm_peak_ids(
                        all_k27ac_peaks,
                        used_peaks,
                        eps,
                        cm_size,
                        self.ncm_k27ac_peak_coord_dict,
                        mark="H3K27ac",
                    )
                else:
                    simulated_cm_peak_ids = self.get_simulated_cm_peak_ids_both_marks(
                        all_k27ac_peaks,
                        all_k4me1_peaks,
                        n_k27ac_peaks,
                        used_peaks,
                        eps,
                        cm_size,
                    )
                if len(simulated_cm_peak_ids) != cm_size:
                    continue

                simulated_cm = sorted(
                    [self.ncm_peak_coord_dict[p_id] for p_id in simulated_cm_peak_ids]
                )
                if not simulated_cm:
                    continue
                if len(self.merge_intersecting_intervals(simulated_cm)) != cm_size:
                    continue
                if len(self.merge_intersecting_intervals(simulated_cm)) == 1:
                    continue

                set_of_all_peaks = set_of_all_peaks - used_peaks
                used_peaks.update(simulated_cm_peak_ids)

                simulated_cm_peak_bp = np.sum(
                    [peak[1] - peak[0] for peak in simulated_cm]
                )
                sim_pid = self.get_pid(simulated_cm)
                sim_mean_mean, _, sim_mean_std, _ = self.get_peak_stats(sim_pid)
                simulated_cms[cm_id + "_" + str(i)] = simulated_cm

                sim_k27ac_fraction = {
                    mark: n_peaks / len(simulated_cm_peak_ids)
                    for mark, n_peaks in Counter(
                        [sim_peak.split(":")[0] for sim_peak in simulated_cm_peak_ids]
                    ).items()
                }.get("H3K27ac", 0.0)
                sim_feature_dict = {
                    "cm_size": cm_size,
                    "cm_length": max(simulated_cm)[1] - min(simulated_cm)[0],
                    "peak_mean": sim_mean_mean,
                    "peak_std": sim_mean_std,
                    "peak_bp_length": simulated_cm_peak_bp,
                    "k27ac_fraction": sim_k27ac_fraction,
                }
                # 'gc_content': self.get_cm_gc_content(sim_pid)}
                sim_feature_vec[cm_id + "_" + str(i)] = sim_feature_dict
                simulated_cms_pids[cm_id + "_" + str(i)] = simulated_cm_peak_ids
            if not simulated_cms:
                continue
            simulated_cms_pid_dict[cm_id] = simulated_cms_pids
            simulated_cms_dict[cm_id] = simulated_cms
            sim_feature_vec_dict[cm_id] = sim_feature_vec

        return (
            ref_feature_vec,
            simulated_cms_dict,
            simulated_cms_pid_dict,
            sim_feature_vec_dict,
        )

    def get_pid(self, peak_coord_list):
        return [
            ":".join([str(self.chromosome)] + list(map(str, p_coord)))
            for p_coord in peak_coord_list
        ]

    def get_peak_stats(self, pid_list):
        if ("H3K" in list(self.chr_count_mtx.index)[0]) or (
            "ATAC" in list(self.chr_count_mtx.index)[0]
        ):
            self.chr_count_mtx.index = self.chr_count_mtx.index.str.replace(
                "^.*?(?=chr)", ""
            )
        pid_in_idx = list(set(self.chr_count_mtx.index).intersection(set(pid_list)))
        peak_stds = (
            self.chr_count_mtx.loc[pid_in_idx, :].dropna(axis=1).values.std(axis=1)
        )
        peak_means = (
            self.chr_count_mtx.loc[pid_in_idx, :].dropna(axis=1).values.mean(axis=1)
        )
        return peak_means.mean(), peak_means.std(), peak_stds.mean(), peak_stds.std()
