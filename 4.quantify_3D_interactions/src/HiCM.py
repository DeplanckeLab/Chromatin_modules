import numpy as np
import pandas as pd

# from numpy import inf


class HiCM:
    def __init__(self, cool_mtx, resolution, chr_sizes_path):
        self.cool_mtx = cool_mtx
        self.resolution = resolution
        self.chr_sizes_path = chr_sizes_path
        chr_sizes_hg19_df = pd.read_csv(self.chr_sizes_path, sep=",")
        self.chr_sizes_hg19 = dict(
            zip(chr_sizes_hg19_df["chr_id"], chr_sizes_hg19_df["chr_length"])
        )

    @staticmethod
    def merge_dags(set1, set2):
        if set1.intersection(set2) != set():
            return set1.union(set2)
        else:
            return None

    @staticmethod
    def prepare_hic(
        numpy_mtx,
        drop_zeros=False,
    ):
        if drop_zeros:
            numpy_mtx = numpy_mtx[~np.all(numpy_mtx == 0, axis=1), :]
            numpy_mtx = numpy_mtx[:, ~np.all(numpy_mtx == 0, axis=0)]
        # if log_norm:
        #     numpy_mtx_log = np.log10(numpy_mtx + 1)
        #     numpy_mtx_log[numpy_mtx_log == -inf] = 0
        #     numpy_mtx = numpy_mtx_log
        return numpy_mtx

    def get_peak_interactions(
        self,
        chromosome,
        vcm_peak_coord,
        drop_zeros=True,
        get_interaction_snap=False,
        n_bins=None,
        window=None,
    ):
        if get_interaction_snap:
            if (n_bins is None) & (window is None):
                window_ = 20000
            elif not n_bins is None:
                window_ = self.resolution * n_bins
            else:
                window_ = window
        else:
            window_ = 0

        vcm_start = min(vcm_peak_coord, key=lambda x: x[0])[0] - window_
        vcm_end = max(vcm_peak_coord, key=lambda x: x[1])[1] + window_

        start = max(int(np.floor(vcm_start / self.resolution) * self.resolution), 0)
        if "chr" in chromosome:
            end = min(
                int(np.ceil(vcm_end / self.resolution) * self.resolution),
                self.chr_sizes_hg19[chromosome] - 1,
            )
        else:
            chromosome = "chr" + chromosome
            end = min(
                int(np.ceil(vcm_end / self.resolution) * self.resolution),
                self.chr_sizes_hg19[chromosome] - 1,
            )
        vcm_idx = [
            [
                np.floor((peak[0] - start) / self.resolution).astype(int),
                np.ceil((peak[1] - start) / self.resolution).astype(int),
            ]
            for peak in vcm_peak_coord
        ]
        if self.cool_mtx.chromnames[:-3]:
            chr_cool = self.cool_mtx.chromnames[:-3][0]
            if "chr" in chr_cool:
                hic_chromosome = chromosome
            else:
                hic_chromosome = chromosome.replace("chr", "")
        hic_mtx = self.cool_mtx.matrix(balance=True, ignore_index=False).fetch(
            hic_chromosome + ":" + str(start) + "-" + str(end)
        )
        final_hic_mtx = self.prepare_hic(hic_mtx, drop_zeros=drop_zeros)
        pairwise_interaction_freq = [
            np.nansum(final_hic_mtx[peak_1[0] : peak_1[1], peak_2[0] : peak_2[1]])
            / final_hic_mtx[peak_1[0] : peak_1[1], peak_2[0] : peak_2[1]].size
            for i, peak_1 in enumerate(vcm_idx)
            for _, peak_2 in enumerate(vcm_idx[i + 1 :])
        ]
        output = (
            np.nansum(pairwise_interaction_freq)
            / np.array(pairwise_interaction_freq).size
        )
        if get_interaction_snap:
            return (output, final_hic_mtx)
        else:
            return output

    def get_vcm_stats(self, chromosome, row, drop_zeros=True):
        start = max(int(np.floor(row["start"] / self.resolution) * self.resolution), 0)
        end = min(
            int(np.ceil(row["end"] / self.resolution) * self.resolution),
            self.chr_sizes_hg19[chromosome] - 1,
        )

        p_starts = (
            np.array([int(p_s) for p_s in row["peak_starts"].split(",")]) + row["start"]
        )
        p_ends = p_starts + np.array(
            [int(p_len) for p_len in row["peak_length"].split(",")]
        )
        nvcm_coord = list(zip(p_starts, p_ends))
        nvcm_idx = [
            [
                np.floor((nDAG_peak[0] - start) / self.resolution).astype(int),
                np.ceil((nDAG_peak[1] - start) / self.resolution).astype(int),
            ]
            for nDAG_peak in nvcm_coord
        ]
        # try:
        hic_mtx = self.cool_mtx.matrix(balance=True, ignore_index=False).fetch(
            chromosome + ":" + str(start) + "-" + str(end)
        )
        #     hic_mtx = hic_mtx[~np.isnan(hic_mtx).all(axis=1)]
        #     hic_mtx = hic_mtx[:, ~np.isnan(hic_mtx).all(axis=0)]

        if drop_zeros:
            hic_mtx_log_no_zeros = self.prepare_hic(hic_mtx, drop_zeros=True)
            np.fill_diagonal(np.array(hic_mtx_log_no_zeros), np.nan)
            final_hic_mtx = hic_mtx_log_no_zeros
        else:
            hic_mtx_log_zeros = self.prepare_hic(hic_mtx, drop_zeros=False)
            np.fill_diagonal(np.array(hic_mtx_log_zeros), np.nan)
            final_hic_mtx = hic_mtx_log_zeros

        pairwise_dd_interaction_freq = [
            np.nanmean(
                final_hic_mtx[
                    nDAG_peak_1[0] : nDAG_peak_1[1], nDAG_peak_2[0] : nDAG_peak_2[1]
                ]
            )
            for i, nDAG_peak_1 in enumerate(nvcm_idx)
            for j, nDAG_peak_2 in enumerate(nvcm_idx[i + 1 :])
        ]
        return np.nanmean(pairwise_dd_interaction_freq)

    def get_dag_stats(self, chromosome, row, dependent_dict, drop_zeros=True):
        start = max(
            int(np.floor(row["start"] / self.resolution) * self.resolution), 0
        )  # - 500000
        end = int(np.ceil(row["end"] / self.resolution) * self.resolution)  # + 500000

        hic_mtx = self.cool_mtx.matrix(balance=True, ignore_index=False).fetch(
            chromosome + ":" + str(start) + "-" + str(end)
        )
        #     hic_mtx = hic_mtx[~np.isnan(hic_mtx).all(axis=1)]
        #     hic_mtx = hic_mtx[:, ~np.isnan(hic_mtx).all(axis=0)]

        if drop_zeros:
            hic_mtx_log_no_zeros = self.prepare_hic(hic_mtx, drop_zeros=True)
            np.fill_diagonal(np.array(hic_mtx_log_no_zeros), np.nan)
            final_hic_mtx = hic_mtx_log_no_zeros
        else:
            hic_mtx_log_zeros = self.prepare_hic(hic_mtx, drop_zeros=False)
            np.fill_diagonal(np.array(hic_mtx_log_zeros), np.nan)
            final_hic_mtx = hic_mtx_log_zeros

        dependent_coord = dependent_dict[row["DAG_id"]]
        dependent_idx = [
            [
                np.floor((dependent[0] - start) / self.resolution).astype(int),
                np.ceil((dependent[1] - start) / self.resolution).astype(int),
            ]
            for dependent in dependent_coord
        ]
        lead_coord = [row["lead_peak_start"], row["lead_peak_end"]]
        lead_idx = [
            np.floor((lead_coord[0] - start) / self.resolution).astype(int),
            np.ceil((lead_coord[1] - start) / self.resolution).astype(int),
        ]

        # all lead-dependent interactions
        pairwise_ld_interaction_freq = [
            np.nanmean(
                final_hic_mtx[lead_idx[0] : lead_idx[1], dependent[0] : dependent[1]]
            )
            for dependent in dependent_idx
        ]
        mean_ld_interaction = np.nanmean(pairwise_ld_interaction_freq)
        if len(dependent_coord) > 1:
            # all dependent-dependent interactions
            pairwise_dd_interaction_freq = [
                np.nanmean(
                    final_hic_mtx[
                        dependent_1[0] : dependent_1[1], dependent_2[0] : dependent_2[1]
                    ]
                )
                for i, dependent_1 in enumerate(dependent_idx)
                for j, dependent_2 in enumerate(dependent_idx[i + 1 :])
            ]
            # all interactions (lead-dependent & dependent-dependent)
            pairwise_ld_interaction_freq.extend(pairwise_dd_interaction_freq)
            mean_dd_interaction = np.nanmean(pairwise_dd_interaction_freq)
        return (
            mean_ld_interaction,
            mean_dd_interaction,
            np.nanmean(pairwise_ld_interaction_freq),
        )

    def get_dags_as_vcms(self, dag_ld_pairs):
        """if method == 'PHM':"""
        dags_union = []
        for i, (lead_peak, peak_set1) in enumerate(dag_ld_pairs):
            temp_set = set(peak_set1)
            lead_set = set([lead_peak])
            if temp_set.issubset(set([el for s in dags_union for el in s[1]])):
                continue
            for j in range(i + 1, len(dag_ld_pairs)):
                if self.merge_dags(temp_set, set(dag_ld_pairs[j][1])) is not None:
                    temp_set.update(set(dag_ld_pairs[j][1]))
                    lead_set.update(set([dag_ld_pairs[j][0]]))
                else:
                    dags_union.append((lead_set, temp_set))
                    break
                continue
        return dags_union


# def get_ep_centered_hic_snap(
#     self,
#     chromosome,
#     vcm_peak_coord,
#     drop_zeros=True,
#     n_bins=4,
#     window=None,
# ):
#     if (n_bins is None) & (window is None):
#         window_ = 20000
#     elif not n_bins is None:
#         window_ = self.resolution * n_bins
#     else:
#         window_ = window
#     vcm_start = min(vcm_peak_coord, key=lambda x: x[0])[0] - window_
#     vcm_end = max(vcm_peak_coord, key=lambda x: x[1])[1] + window_
#     start = max(int(np.floor(vcm_start / self.resolution) * self.resolution), 0)
#     if "chr" in chromosome:
#         end = min(
#             int(np.ceil(vcm_end / self.resolution) * self.resolution),
#             self.chr_sizes_hg19[chromosome] - 1,
#         )
#     else:
#         end = min(
#             int(np.ceil(vcm_end / self.resolution) * self.resolution),
#             self.chr_sizes_hg19["chr" + chromosome] - 1,
#         )
#     hic_mtx = self.cool_mtx.matrix(balance=True, ignore_index=False).fetch(
#         "chr" + chromosome + ":" + str(start) + "-" + str(end)
#     )
#     return self.prepare_hic(hic_mtx, drop_zeros=drop_zeros)
