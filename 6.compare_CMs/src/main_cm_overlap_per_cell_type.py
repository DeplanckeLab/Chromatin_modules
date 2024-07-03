import argparse
import itertools
import logging
import numpy as np
import os
import pandas as pd

from CM_intersection_class import CM_Overlap
from utils import handle_threshold, process_thresholds


parser = argparse.ArgumentParser(
    description="Chromatin module similarity scores ©Olga Pushkarev"
)
parser.add_argument(
    "-m", "--method", type=str, help="Input method (or dataset)", required=True
)
parser.add_argument(
    "-d",
    "--dataset_list",
    type=str,
    help="Comma-separated input list of datasets (or methods)",
    required=True,
)
parser.add_argument(
    "-cm", "--cm_path", type=str, help="Path to chromatin modules", required=True
)
parser.add_argument(
    "-pcm",
    "--peak_path",
    type=str,
    help="Path to chromatin module peaks",
    required=True,
)
# phm parameters
parser.add_argument(
    "-pp",
    "--pp_threshold",
    help="Posterior probability of causal intercations for phm",
    default="0.8",
    required=True,
)
# VCMtools parameters
parser.add_argument(
    "-vw",
    "--vcm_window",
    help="cis- window used in VCMtools",
    default="0.5Mb",
    type=str,
    required=True,
)
parser.add_argument(
    "-pv",
    "--pv_threshold",
    help="p-value threshold used in VCMtools",
    default=0.001,
    required=True,
)
# Clomics parameters
parser.add_argument(
    "-np",
    "--n_peaks",
    help="Number of peaks used in Clomics",
    default="n_peaks_200",
    type=str,
    required=True,
)
parser.add_argument(
    "-bg",
    "--bg_threshold",
    help="Background correlation threshold used in Clomics",
    default="bg_threshold_3",
    type=str,
    required=True,
)
# Peak overlap path
parser.add_argument(
    "-pop",
    "--peak_overlap_path",
    type=str,
    help="Path to overlapping peaks",
    required=True,
)
parser.add_argument(
    "-q",
    "--use_cmQTLs",
    type=int,
    help="Use only chromatin modules with QTLs",
    default=0,
    required=False,
)
parser.add_argument(
    "-mcs",
    "--min_cm_size",
    type=int,
    help="Use only chromatin modules of size >= min_cm_size",
    default=0,
    required=False,
)
parser.add_argument(
    "-s",
    "--score_type",
    type=str,
    default="HM",
    help="Average type: either HM (harmonic mean) or AM (Arithmetic mean)",
    required=True,
)
parser.add_argument(
    "-t", "--similarity_threshold", help="Description for bar argument", required=True
)
parser.add_argument(
    "-o",
    "--output_path",
    type=str,
    help="Output directory for the scores",
    required=True,
)
parser.add_argument(
    "-ow",
    "--overwrite",
    type=int,
    help="Output directory for the scores",
    default=1,
    required=False,
)
args = parser.parse_args()

# Arguments
dataset = args.method
methods_list = args.dataset_list.split(",")
core_cm_path = args.cm_path
pp_threshold = args.pp_threshold
vcm_window = args.vcm_window
pv_threshold = args.pv_threshold
n_peaks = args.n_peaks
bg_threshold = args.bg_threshold

path_to_input_peak_directory = args.peak_path
path_to_input_peak_overlap_directory = args.peak_overlap_path
use_only_vcms_with_QTL = args.use_cmQTLs
min_cm_size = int(args.min_cm_size)
score_type = args.score_type
similarity_threshold = args.similarity_threshold
path_to_output_directory = args.output_path
overwrite = args.overwrite

logging.basicConfig(
    filename=os.path.join(path_to_output_directory, "logfile.txt"),
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    filemode="a",
)
logging.info("######################################################")
logging.info("# Chromatin module similarity scores ©Olga Pushkarev #")
logging.info("######################################################")
logging.info("\nList of arguments:\n")
logging.info("Input:")
logging.info("\t 1. Input dataset: " + dataset)
logging.info("\t 2. Using the following methods: " + ", ".join(methods_list))
logging.info("\t 3. Path to chromatin module peaks: " + path_to_input_peak_directory)
logging.info("\t 4. Path to overlapping peaks: " + path_to_input_peak_overlap_directory)
if use_only_vcms_with_QTL:
    cmQTLs_bool = "Yes"
else:
    cmQTLs_bool = "No"
logging.info("\t 5. Use only chromatin modules with cmQTLs? " + cmQTLs_bool)
if score_type == "HM":
    score_type_name = "HM ( Harmonic mean )"
else:
    score_type_name = "AM ( Arithmetic mean )"
logging.info("\t 6. Score type: " + score_type_name)
if "-" in similarity_threshold:
    range_start = float(similarity_threshold.split("-")[0])
    range_end = float(similarity_threshold.split("-")[1].split(";")[0])
    step = float(similarity_threshold.split("-")[1].split(";")[1])
    logging.info(
        "".join(
            [
                "\t 7. Similarity thresholds in range [",
                str(range_start),
                ", ",
                str(range_end),
                "] with step ",
                str(step),
            ]
        )
    )
else:
    logging.info("\t 7. Similarity threshold(s): " + similarity_threshold)
if overwrite:
    overwrite_bool = "Yes"
else:
    overwrite_bool = "No"
logging.info("\t 8. Overwrite .npy files with scores?: " + overwrite_bool + "\n")
logging.info("\t 9. Posterior probability threshold for phm: " + pp_threshold + "\n")
logging.info("Output:")
logging.info(
    "\t 1. Output directory for similarity scores:" + path_to_output_directory + "\n"
)

path_to_input_peak_directory = path_to_input_peak_directory.replace("\\", "/")
path_to_input_peak_overlap_directory = path_to_input_peak_overlap_directory.replace(
    "\\", "/"
)
path_to_output_directory = os.path.join(
    path_to_output_directory.replace("\\", "/")
    # 'pp_threshold_' + pp_threshold
)
if not os.path.exists(path_to_output_directory):
    os.makedirs(path_to_output_directory)

crd_tracks_path = os.path.join(
    core_cm_path,
    "clomics",
    n_peaks,
    bg_threshold,
    dataset,
    dataset + "_Clomics_CM.tracks.bed",
)
vcm_tracks_path = os.path.join(
    core_cm_path,
    "vcmtools",
    "VCMs",
    vcm_window,
    dataset,
    "_".join(
        [
            dataset
            + "_NOT_TOTEM_VCMs_corrected_pvalue_"
            + str(pv_threshold)
            + ".tracks.bed"
        ]
    ),
)
phm_tracks_path = os.path.join(
    core_cm_path,
    "phm",
    dataset + "_phm_tracks_content",
    "_".join([dataset, pp_threshold, "merged_phm_all_chr_NOT_TOTEM.tracks.bed"]),
)

cm_overlap = CM_Overlap(
    output_path=path_to_output_directory,
    methods_list=methods_list,
    score_type=score_type,
    use_cmQTLs=use_only_vcms_with_QTL,
)

crds_with_qtl = None
vcms_with_qtl = None
if use_only_vcms_with_QTL:
    crds_with_qtl = cm_overlap.get_cms_with_cmQTL(core_cm_path, "clomics", dataset)
    vcms_with_qtl = cm_overlap.get_cms_with_cmQTL(core_cm_path, "vcmtools", dataset)

crd_peaks_df = cm_overlap.get_cm_peaks(
    os.path.join(path_to_input_peak_directory, dataset + "_clomics_all_peaks.bed")
)
vcm_peaks_df = cm_overlap.get_cm_peaks(
    os.path.join(
        path_to_input_peak_directory,
        dataset + "_vcmtools_NOT_TOTEM_peaks.bed",
    )
)
phm_peaks_df = cm_overlap.get_cm_peaks(
    os.path.join(
        path_to_input_peak_directory,
        dataset + "_phm_NOT_TOTEM_peaks.bed",
    )
)

final_crds_with_qtl = None
final_vcms_with_qtl = None
if min_cm_size != 0:
    crd_ids_min_size = cm_overlap.subset_cms_by_size(crd_tracks_path, min_cm_size)
    vcm_ids_min_size = cm_overlap.subset_cms_by_size(vcm_tracks_path, min_cm_size)
    phm_ids_min_size = cm_overlap.subset_cms_by_size(phm_tracks_path, min_cm_size)

    # crd_ids_min_size subset tracks and peaks
    final_crd_peaks_df = crd_peaks_df[crd_peaks_df.iloc[:, 4].isin(crd_ids_min_size)]
    final_vcm_peaks_df = vcm_peaks_df[vcm_peaks_df.iloc[:, 4].isin(vcm_ids_min_size)]
    final_phm_peaks_df = phm_peaks_df[phm_peaks_df.iloc[:, 4].isin(phm_ids_min_size)]

    min_sizes_dict = {
        "clomics": crd_ids_min_size,
        "vcmtools": vcm_ids_min_size,
        "phm": phm_ids_min_size,
    }

    if use_only_vcms_with_QTL:
        final_crds_with_qtl = list(set(crds_with_qtl).intersection(crd_ids_min_size))
        final_vcms_with_qtl = list(set(vcms_with_qtl).intersection(vcm_ids_min_size))
else:
    final_crd_peaks_df = crd_peaks_df
    final_vcm_peaks_df = vcm_peaks_df
    final_phm_peaks_df = phm_peaks_df

    if use_only_vcms_with_QTL:
        final_crds_with_qtl = crds_with_qtl
        final_vcms_with_qtl = vcms_with_qtl

del crd_peaks_df, vcm_peaks_df, phm_peaks_df
del crds_with_qtl, vcms_with_qtl

method_annotation_dict = {
    "clomics": {
        "peaks_df": final_crd_peaks_df,
        "cm_qtl_list": final_crds_with_qtl,
        "peak_to_cm_dict": dict(
            zip(
                list(final_crd_peaks_df.iloc[:, 3]), list(final_crd_peaks_df.iloc[:, 4])
            )
        ),
    },
    "vcmtools": {
        "peaks_df": final_vcm_peaks_df,
        "cm_qtl_list": final_vcms_with_qtl,
        "peak_to_cm_dict": dict(
            zip(
                list(final_vcm_peaks_df.iloc[:, 3]), list(final_vcm_peaks_df.iloc[:, 4])
            )
        ),
    },
    "phm": {
        "peaks_df": final_phm_peaks_df,
        "cm_qtl_list": None,
        "peak_to_cm_dict": dict(
            zip(
                list(final_phm_peaks_df.iloc[:, 3]), list(final_phm_peaks_df.iloc[:, 4])
            )
        ),
    },
}
if overwrite:
    logging.info("Calculating the similarity scores for chromatin modules...")
    all_method_pairs = list(itertools.permutations(methods_list, 2)) + [
        [method, method] for method in methods_list
    ]
    for method_A, method_B in all_method_pairs:
        A_B_peak_overlap = pd.read_csv(
            os.path.join(
                path_to_input_peak_overlap_directory,
                "_".join([method_A, method_B, "jaccard_peak_overlap.bed"]),
            ),
            sep="\t",
            header=None,
        )
        A_B_peak_overlap.iloc[:, 4] = A_B_peak_overlap.iloc[:, 4].str.replace("_", "~")
        A_B_peak_overlap.iloc[:, 10] = A_B_peak_overlap.iloc[:, 10].str.replace(
            "_", "~"
        )
        if min_cm_size != 0:
            A_B_peak_overlap = A_B_peak_overlap[
                (A_B_peak_overlap.iloc[:, 4].isin(min_sizes_dict[method_A]))
                & (A_B_peak_overlap.iloc[:, 10].isin(min_sizes_dict[method_B]))
            ]
        (
            A_B_peak_overlap_df,
            A_B_peak_overlap_dict,
        ) = cm_overlap.get_jaccard_peak_overlap(
            A_B_peak_overlap,
            method_A_peak_start_col=1,
            method_A_peak_end_col=2,
            method_A_vcm_id_col=4,
            method_B_peak_start_col=7,
            method_B_peak_end_col=8,
            method_B_vcm_id_col=10,
            peak_len_overlap_col=12,
            A_CMs_with_qtl=method_annotation_dict[method_A]["cm_qtl_list"],
            B_CMs_with_qtl=method_annotation_dict[method_B]["cm_qtl_list"],
            bp_jaccard=False,
        )
        A_B_F1_scores = cm_overlap.get_scores_for_methods_AB_dict(
            dataset=dataset,
            method_A=method_A,
            method_B=method_B,
            method_A_peaks_df=method_annotation_dict[method_A]["peaks_df"],
            method_B_peaks_df=method_annotation_dict[method_B]["peaks_df"],
            method_AB_peak_overlap_df=A_B_peak_overlap_df,
            method_AB_overlapping_modules_dict=A_B_peak_overlap_dict,
            A_CMs_with_qtl=method_annotation_dict[method_A]["cm_qtl_list"],
            B_CMs_with_qtl=method_annotation_dict[method_B]["cm_qtl_list"],
            method_A_cm_id_col_in_peaks=4,
            method_A_peak_id_col_in_overlap=3,
            method_A_cm_id_col_in_overlap=4,
            method_B_cm_id_col_in_peaks=4,
            method_B_peak_id_col_in_overlap=9,
            method_B_cm_id_col_in_overlap=10,
            save_scores=True,
        )

similarity_thresholds = handle_threshold(similarity_threshold)

for threshold in similarity_thresholds:
    threshold = np.round(threshold, 3)
    if threshold == 0:
        direction_list = ["greater"]
    else:
        direction_list = ["less", "greater"]

    universal_crds_tracks_le, universal_crds_tracks_gr = process_thresholds(
        cm_overlap, "clomics", dataset, crd_tracks_path, threshold, min_cm_size
    )
    universal_vcms_tracks_le, universal_vcms_tracks_gr = process_thresholds(
        cm_overlap, "vcmtools", dataset, vcm_tracks_path, threshold, min_cm_size
    )
    universal_phms_tracks_le, universal_phms_tracks_gr = process_thresholds(
        cm_overlap, "phm", dataset, phm_tracks_path, threshold, min_cm_size
    )
