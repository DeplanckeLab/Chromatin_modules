import argparse
import logging
import numpy as np
import os
import pandas as pd

from DAG2CM_class import DAG2CM

parser = argparse.ArgumentParser(
    description="Build CM-like track and content files from DAGs ©Olga Pushkarev"
)
parser.add_argument("-d", "--dataset", type=str, help="Input dataset", required=True)
parser.add_argument(
    "-i",
    "--path_to_input_directory",
    type=str,
    help="Path to folder with count matrices",
    default="path_to_output",
    required=False,
)
parser.add_argument(
    "-o",
    "--path_to_output_directory",
    type=str,
    help="Path to output folder",
    default=0,
    required=True,
)
parser.add_argument(
    "-t",
    "--pp_threshold",
    type=str,
    help="Comma separated posterior probability thresholds for causal peak interaction hypothesis, e.g., 0.5,0.6,0.7 or 0.5",
    default="0.5",
    required=False,
)
parser.add_argument(
    "-c", "--chromosomes", type=str, help="Chromosomes to process", required=True
)

args = parser.parse_args()

# Arguments
dataset = args.dataset
pp_threshold = float(args.pp_threshold)
chromosomes_str = args.chromosomes

# +
path_to_input_directory = args.path_to_input_directory
path_to_output_directory = args.path_to_output_directory

if path_to_input_directory == "path_to_output":
    path_to_input_directory = path_to_output_directory
# -

if "," in chromosomes_str:
    if not "chr" in chromosomes_str.split(",")[0]:
        chromosomes = ["chr" + str(chr_id) for chr_id in chromosomes_str.split(",")]
    else:
        chromosomes = chromosomes_str.split(",")
elif "-" in chromosomes_str:
    chromosomes_lst = [el.replace("chr", "") for el in chromosomes_str.split("-")]
    chromosomes = [
        "chr" + str(chr_int)
        for chr_int in np.arange(int(chromosomes_lst[0]), int(chromosomes_lst[1]) + 1)
    ]
else:
    if not "chr" in chromosomes_str:
        chromosomes = ["chr" + chromosomes_str]
    else:
        chromosomes = [chromosomes_str]

logging.basicConfig(
    filename=os.path.join(path_to_output_directory, dataset + "_logfile.txt"),
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    filemode="a",
)
logging.info("####################################################################")
logging.info("# Build CM-like track and content files from DAGs ©Olga Pushkarev #")
logging.info("####################################################################")
logging.info("\nList of arguments:\n")
logging.info("Input:")
logging.info("\t 1. Input dataset: " + dataset)
path_to_input_directory = path_to_input_directory.replace("\\", "/")
path_to_output_directory = path_to_output_directory.replace("\\", "/")
logging.info("\t 2. Path to input files): " + path_to_input_directory)
logging.info("\n")
logging.info("Output:")
logging.info("\t 1. Path to output files: " + path_to_output_directory)
logging.info("\n")

# # Set working directories
tracks_content_path = os.path.join(
    path_to_output_directory, "_".join([dataset, "phm_tracks_content"])
)
DAGs_path = os.path.join(path_to_output_directory, "_".join([dataset, "DAGs"]))
if not os.path.exists(tracks_content_path):
    os.makedirs(tracks_content_path)
if not os.path.exists(DAGs_path):
    os.makedirs(DAGs_path)
DAG2CM_class = DAG2CM(path_to_output_directory, chromosomes, pp_threshold, dataset)
logging.info("Merging data from all chromosomes...")
merged_tracks = DAG2CM_class.merge_tracks()
merged_contents = DAG2CM_class.merge_contents()
merged_dags = DAG2CM_class.merge_dags()

logging.info("Saving merged files...")
merged_tracks.to_csv(
    os.path.join(
        tracks_content_path,
        "_".join([dataset, str(pp_threshold), "merged_phm_all_chr.tracks.bed"]),
    ),
    sep="\t",
    index=False,
    header=False,
)
merged_contents.to_csv(
    os.path.join(
        tracks_content_path,
        "_".join([dataset, str(pp_threshold), "merged_phm_all_chr.content.txt"]),
    ),
    sep="\t",
    index=False,
    header=False,
)
merged_dags.to_csv(
    os.path.join(
        DAGs_path, "_".join([str(pp_threshold), "merged_DAGs_for_all_chr.tsv"])
    ),
    sep="\t",
    index=False,
    header=False,
)
if merged_tracks.shape[0] != 0:
    tracks_column_names = [
        "chr",
        "start",
        "end",
        "cm_id",
        "number",
        "strain",
        "start_duplicate",
        "end_duplicate",
        "numbers",
        "cm_size",
        "peak_length",
        "peak_starts",
        "is_totem",
    ]
    merged_tracks.columns = tracks_column_names
    merged_tracks[merged_tracks["is_totem"] == 0].to_csv(
        os.path.join(
            tracks_content_path,
            "_".join(
                [dataset, str(pp_threshold), "merged_phm_all_chr_NOT_TOTEM.tracks.bed"]
            ),
        ),
        sep="\t",
        index=False,
        header=False,
    )
    contents_column_names = ["cm_id", "n_clusters", "peaks"]
    merged_contents.columns = contents_column_names
    not_totem_contents = merged_contents[
        merged_contents["cm_id"].isin(
            merged_tracks[merged_tracks["is_totem"] == 0]["cm_id"]
        )
    ]
    not_totem_contents.to_csv(
        os.path.join(
            tracks_content_path,
            "_".join(
                [dataset, str(pp_threshold), "merged_phm_all_chr_NOT_TOTEM.content.txt"]
            ),
        ),
        sep="\t",
        index=False,
        header=False,
    )

# Save files that do not dependent on pp_threshold
peak_coord_path = os.path.join(
    path_to_output_directory, "_".join([dataset, "all_chr_peak_coordinates.tsv"])
)
lead_variants_path = os.path.join(
    path_to_output_directory,
    "_".join([dataset, "ALL_lead_caQTL_variants_merged_chr.tsv"]),
)

if not os.path.exists(peak_coord_path):
    logging.info("File with peaks is missing! Preparing and saving peak coordinates...")
    # Merge peak files
    peaks_to_merge = []
    for chromosome in chromosomes:
        peaks_chr = pd.read_csv(
            os.path.join(
                path_to_input_directory,
                chromosome,
                "_".join([chromosome, "peak_coordinates.bed.gz"]),
            ),
            sep="\t",
            header=None,
        ).reset_index(drop=True)
        peaks_chr["peak_id"] = np.arange(1, peaks_chr.shape[0] + 1)
        peaks_chr["peak_id"] = peaks_chr["peak_id"].astype(str) + "_" + chromosome
        peaks_to_merge.append(peaks_chr)
    merged_peaks = pd.concat(peaks_to_merge, axis=0)
    merged_peaks.to_csv(peak_coord_path, sep="\t", index=False, header=False)

# I would recommend to start with finding the naive lead variant by merging the input file and the posterior
# probability “pp.gz" given by hm. The 4th column of pp.gz gives the posterior probability of a variant being causal
# QTL for each peak. You can simply pick up one gives the highest value within each peak
# (stratified by the first column of pp.gz).

# +
# if not os.path.exists(lead_variants_path):
#     # Select lead variants
#     logging.info(
#         "File with lead variants is missing! Preparing and saving lead variants..."
#     )
#     for chromosome in chromosomes:
#         all_variants_df = pd.read_csv(
#             os.path.join(path_to_output_directory, chromosome, "bayeslm_1.gz"),
#             sep="\t",
#             header=None,
#         )
#         all_variants_df.columns = [
#             "Peak",
#             "Chr",
#             "Pos",
#             "RsID",
#             "Ref",
#             "Alt",
#             "AF",
#             "R2",
#             "Variant_loc",
#             "Inside_Peak",
#             "Beta",
#             "SE",
#             "Log10_BF",
#         ]
#         hm_pp = pd.read_csv(
#             os.path.join(path_to_output_directory, chromosome, "hm_output", "pp.gz"),
#             sep="\t",
#             header=None,
#         )
#         hm_pp.columns = ["Peak", "X", "XX", "pp_of_causal_QTL"]
#         hm_pp = hm_pp[["Peak", "pp_of_causal_QTL"]]

#         peaks_and_lead_snps = pd.DataFrame(
#             columns=[
#                 "Peak",
#                 "chr",
#                 "rs_pos",
#                 "RsID",
#                 "Ref_Alt",
#                 "AF",
#                 "R2",
#                 "rs_location",
#                 "inside_peak",
#                 "Beta",
#                 "SE",
#                 "Log10_BF",
#                 "pp_of_caQTL",
#             ],
#             index=np.arange(len(set(all_variants_df["Peak"]))),
#         )
#         for idx, (peak_id, peak_df) in enumerate(all_variants_df.groupby("Peak")):
#             hm_pp_chr = hm_pp[hm_pp["Peak"] == peak_id]
#             if hm_pp_chr.shape[0] != peak_df.shape[0]:
#                 continue
#                 # logging.info('ERROR: Shapes of peak and hm_pp dataframes are different!')
#                 # sys.exit('ERROR! Check shapes of peak and hm_pp dataframes')
#                 # break
#             rs_peaks_chr = pd.concat([peak_df, hm_pp_chr], axis=1)

#             max_pp_variants = rs_peaks_chr[
#                 rs_peaks_chr["pp_of_causal_QTL"]
#                 == rs_peaks_chr["pp_of_causal_QTL"].max()
#             ]
#             if max_pp_variants.shape[0] == 0:
#                 continue
#             rs_coord_pair = sorted(
#                 list(zip(max_pp_variants["RsID"].astype(str), max_pp_variants["Pos"]))
#             )
#             rs_ids = ";".join([rs[0] for rs in rs_coord_pair])
#             coordinates = ";".join([str(pos[1]) for pos in rs_coord_pair])

#             rsloc_inpeak_pair = sorted(
#                 list(
#                     zip(max_pp_variants["Variant_loc"], max_pp_variants["Inside_Peak"])
#                 )
#             )
#             rs_loc = ";".join([str(rloc[0]) for rloc in rsloc_inpeak_pair])
#             in_peak = ";".join([str(ip[1]) for ip in rsloc_inpeak_pair])

#             ref_alt_pair = sorted(
#                 list(zip(max_pp_variants["Ref"], max_pp_variants["Alt"]))
#             )
#             ref_alt_str = ";".join(
#                 [str(ref_alt[0]) + ":" + str(ref_alt[1]) for ref_alt in ref_alt_pair]
#             )

#             peaks_and_lead_snps["Peak"].loc[idx] = (
#                 str(peak_id) + "_chr" + str(list(max_pp_variants["Chr"].values)[0])
#             )
#             peaks_and_lead_snps["chr"].loc[idx] = list(max_pp_variants["Chr"].values)[0]
#             peaks_and_lead_snps["rs_pos"].loc[idx] = coordinates
#             peaks_and_lead_snps["RsID"].loc[idx] = rs_ids
#             peaks_and_lead_snps["Ref_Alt"].loc[idx] = ref_alt_str
#             peaks_and_lead_snps["AF"].loc[idx] = list(max_pp_variants["AF"].values)[0]
#             peaks_and_lead_snps["R2"].loc[idx] = list(max_pp_variants["R2"].values)[0]
#             peaks_and_lead_snps["rs_location"].loc[idx] = rs_loc
#             peaks_and_lead_snps["inside_peak"].loc[idx] = in_peak
#             peaks_and_lead_snps["Beta"].loc[idx] = list(max_pp_variants["Beta"].values)[
#                 0
#             ]
#             peaks_and_lead_snps["SE"].loc[idx] = list(max_pp_variants["SE"].values)[0]
#             peaks_and_lead_snps["Log10_BF"].loc[idx] = list(
#                 max_pp_variants["Log10_BF"].values
#             )[0]
#             peaks_and_lead_snps["pp_of_caQTL"].loc[idx] = list(
#                 max_pp_variants["pp_of_causal_QTL"].values
#             )[0]

#         peaks_and_lead_snps.to_csv(
#             os.path.join(
#                 path_to_output_directory,
#                 chromosome,
#                 "_".join([dataset, "ALL_lead_caQTL_variants.tsv"]),
#             ),
#             sep="\t",
#         )
#         lead_variants_df = pd.read_csv(
#             os.path.join(path_to_output_directory, chromosome, "bayeslm_1.gz"),
#             sep="\t",
#             header=None,
#         )
#         lead_variants_df.columns = [
#             "Peak",
#             "Chr",
#             "Pos",
#             "RsID",
#             "Ref",
#             "Alt",
#             "AF",
#             "R2",
#             "P_Lead",
#             "Inside_Peak",
#             "Beta",
#             "SE",
#             "Log10_BF",
#         ]
#         lead_variants_unique = lead_variants_df.loc[
#             lead_variants_df.groupby("Peak")["Log10_BF"].idxmax()
#         ]
#         del lead_variants_df
#         hm_pp = pd.read_csv(
#             os.path.join(path_to_output_directory, chromosome, "hm_output", "pp.gz"),
#             sep="\t",
#             header=None,
#         )
#         hm_pp.columns = ["peak_id", "X", "XX", "pp_of_causal_QTL"]
#         hm_pp_max_pp = (
#             hm_pp.groupby("peak_id")
#             .agg({"pp_of_causal_QTL": "max"})["pp_of_causal_QTL"]
#             .reset_index()
#         )
#         lead_variants_unique.reset_index(drop=True, inplace=True)
#         hm_pp_max_pp.reset_index(drop=True, inplace=True)
#         variant_peaks = pd.concat([lead_variants_unique, hm_pp_max_pp], axis=1)
#         variant_peaks.rename(columns={"P_Lead": "Variant_loc"}, inplace=True)
#         del lead_variants_unique
#         del hm_pp_max_pp
#         variant_peaks.rename(columns={"P_Lead": "Variant_loc"}, inplace=True)
#         variant_peaks["Peak"] = variant_peaks["Peak"].astype(str) + "_" + chromosome
#         variant_peaks.to_csv(
#             os.path.join(
#                 path_to_output_directory,
#                 chromosome,
#                 "_".join([dataset, "lead_caQTL_variants.tsv"]),
#             ),
#             sep="\t",
#         )
#     dfs = [
#         pd.read_csv(
#             os.path.join(
#                 path_to_output_directory,
#                 chromosome,
#                 "_".join([dataset, "ALL_lead_caQTL_variants.tsv"]),
#             ),
#             sep="\t",
#             index_col=0,
#         )
#         for chromosome in chromosomes
#     ]
#     pd.concat(dfs).to_csv(
#         os.path.join(
#             path_to_output_directory,
#             "_".join([dataset, "ALL_lead_caQTL_variants_merged_chr.tsv"]),
#         ),
#         sep="\t",
#         index=False,
#     )
