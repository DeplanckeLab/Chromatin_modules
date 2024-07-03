# %%
import numpy as np
import os
import pandas as pd
import scipy.stats as stats
import matplotlib
import matplotlib.pyplot as plt
import tqdm.auto as tqdm

import collections
import rpy2
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()

# %%
save_files = True
core_path = os.getcwd()

methods = ["vcmtools", "clomics", "phm"]
dataset = "test_data"

peak_path = "/data/pushkare/Chromatin_modules/2.peaks_in_CMs/peak_files"
color_dict = {"H3K4me1": "#afafd8", "H3K27ac": "#cdd6ae"}

# %% [markdown]
# ### Data preparation

# %%
all_peaks_df = pd.read_csv(
    "/data/pushkare/Chromatin_modules/0.input_peaks/peak_bed_files/k27ac_k4me1_peaks.bed",
    sep="\t",
    header=None,
)
all_peaks_df.columns = ["seqnames", "start", "end", "peak_id", "strand"]

# %% [markdown]
# ### Extract sequences for peaks and calculate GC content for them

# %%
bs_genome = importr("BSgenome.Hsapiens.UCSC.hg19")
bio_strings = importr("Biostrings")
genomic_ranges = importr("GenomicRanges")

all_peaks_df = all_peaks_df.set_index("peak_id")
# %%
rpy2.robjects.globalenv["all_peaks_df"] = pandas2ri.py2rpy(all_peaks_df)
granges_for_peaks = genomic_ranges.makeGRangesFromDataFrame(
    rpy2.robjects.globalenv["all_peaks_df"]
)
sequences_for_peaks = bio_strings.getSeq(bs_genome.Hsapiens, granges_for_peaks)
sequences_for_peaks_char = rpy2.robjects.r["as.character"](sequences_for_peaks)
r_dict_with_sequences_for_peaks = dict(
    zip(all_peaks_df.index, sequences_for_peaks_char)
)
# %%
GC_peak_dict = {}
for peak_id, sequence in tqdm.tqdm(r_dict_with_sequences_for_peaks.items()):
    CG_sequence_dict = dict(collections.Counter(sequence))
    if not "G" in sequence:
        CG_sequence_dict["G"] = 0
    if not "C" in sequence:
        CG_sequence_dict["C"] = 0
    GC_peak_dict[peak_id] = (CG_sequence_dict["C"] + CG_sequence_dict["G"]) / len(
        sequence
    )
print("DONE")

# %%
if save_files:
    np.save(
        os.path.join(core_path, "k27ac_k4me1_peaks_GC_content.npy"),
        GC_peak_dict,
        allow_pickle=True,
    )

# %%
GC_peak_dict = np.load(
    os.path.join(core_path, "k27ac_k4me1_peaks_GC_content.npy"),
    allow_pickle=True,
).item()


## Get average GC content per chromatin module
# %%
list_with_gc_per_method = []
dict_with_mean_gc_per_method = {}
for i, method in enumerate(methods):
    if method != "clomics":
        not_totem_cm_peaks = (
            pd.read_csv(
                os.path.join(
                    peak_path, "_".join([dataset, method, "NOT_TOTEM_peaks.bed"])
                ),
                sep="\t",
                header=None,
            )
            .iloc[:, 3]
            .to_list()
        )
        not_totem_peaks_gc_content_df = pd.DataFrame.from_dict(
            {nt_p_id: GC_peak_dict.get(nt_p_id) for nt_p_id in not_totem_cm_peaks},
            orient="index",
            columns=[dataset],
        )
    all_cm_peaks = (
        pd.read_csv(
            os.path.join(peak_path, "_".join([dataset, method, "all_peaks.bed"])),
            sep="\t",
            header=None,
        )
        .iloc[:, 3]
        .to_list()
    )
    not_cm_peaks = (
        pd.read_csv(
            os.path.join(peak_path, "_".join([dataset, "not", method, "peaks.bed"])),
            sep="\t",
            header=None,
        )
        .iloc[:, 3]
        .to_list()
    )

    all_peaks_gc_content_df = pd.DataFrame.from_dict(
        {all_p_id: GC_peak_dict.get(all_p_id) for all_p_id in all_cm_peaks},
        orient="index",
        columns=[dataset],
    )
    not_cm_peaks_gc_content_df = pd.DataFrame.from_dict(
        {not_p_id: GC_peak_dict.get(not_p_id) for not_p_id in not_cm_peaks},
        orient="index",
        columns=[dataset],
    )
    if save_files:
        if not os.path.exists(os.path.join(core_path, "output", method)):
            os.makedirs(os.path.join(core_path, "output", method))
        all_peaks_gc_content_df.to_csv(
            os.path.join(
                core_path,
                "output",
                method,
                "_".join(
                    [method, "GC_content_in", dataset, "per_peak_for_all_peaks.txt"]
                ),
            ),
            sep="\t",
            index=False,
            header=False,
        )
        not_cm_peaks_gc_content_df.to_csv(
            os.path.join(
                core_path,
                "output",
                method,
                "_".join(
                    [
                        method,
                        "GC_content_in",
                        dataset,
                        "per_peak_for_not",
                        method,
                        "peaks.txt",
                    ]
                ),
            ),
            sep="\t",
            index=False,
            header=False,
        )
        if method != "clomics":
            not_totem_peaks_gc_content_df.to_csv(
                os.path.join(
                    core_path,
                    "output",
                    method,
                    "_".join(
                        [
                            method,
                            "GC_content_in",
                            dataset,
                            "per_peak_for_NOT_TOTEM_peaks_all_datasets.txt",
                        ]
                    ),
                ),
                sep="\t",
                index=False,
                header=False,
            )
    dict_with_mean_gc_per_method[method] = {
        "cm_peaks": np.array(list(all_peaks_gc_content_df[dataset].values)) * 100,
        "not_cm_peaks": np.array(list(not_cm_peaks_gc_content_df[dataset].values))
        * 100,
    }
    list_with_gc_per_method.append(
        [
            method,
            np.mean(all_peaks_gc_content_df[dataset].values),
            np.mean(not_cm_peaks_gc_content_df[dataset].values),
        ]
    )
mean_gc_per_ds = [
    dict_with_mean_gc_per_method,
    pd.DataFrame(
        list_with_gc_per_method,
        columns=["method", "GC content\nCM peaks", "GC content\nnon-CM peaks"],
    ).set_index("method")
    * 100,
]

# %%
gc_per_method, _ = mean_gc_per_ds
signif_diff_dict_per_method = {}
for j, (method, gc_dict) in enumerate(gc_per_method.items()):
    cm_peaks = gc_dict["cm_peaks"]
    not_cm_peaks = gc_dict["not_cm_peaks"]

    dict_with_signif_pv = {}
    _, cm_not_cm_pv = stats.mannwhitneyu(cm_peaks, not_cm_peaks)
    if cm_not_cm_pv <= 0.05:
        dict_with_signif_pv["cm-not_cm"] = cm_not_cm_pv

    signif_diff_dict_per_method[method] = dict_with_signif_pv
