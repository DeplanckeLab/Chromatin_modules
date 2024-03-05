import argparse
import logging
import os
import pandas as pd

parser = argparse.ArgumentParser(description="Convert CRDs to CMs ©Olga Pushkarev")
parser.add_argument(
    "-d", "--dataset_name", type=str, help="Input dataset", required=True
)
parser.add_argument(
    "-o",
    "--path_to_cm_batch_files",
    type=str,
    help="Output directory for CM tracks and content files",
    required=True,
)

args = parser.parse_args()

# Arguments
dataset_name = args.dataset_name
path_to_output_directory = args.path_to_cm_batch_files
path_to_output_directory = path_to_output_directory.replace("\\", "/")

logging.basicConfig(
    filename=os.path.join(path_to_output_directory, dataset_name + "_logfile.txt"),
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    filemode="a",
)
logging.info("########################################")
logging.info("# Merge CM batches to one batch ©Olga Pushkarev #")
logging.info("########################################")
logging.info("\nList of arguments:\n")
logging.info("Input:")
logging.info("\t 1. Input dataset: " + dataset_name)
logging.info("\nOutput:")
logging.info(
    "\t 1. Output directory for merged and formated batches of CRDs in CM-like format, i.e. tracks and content files: "
    + path_to_output_directory
)
logging.info("\n")

path = os.path.join(
    path_to_output_directory,
    "_".join([dataset_name, "Clomics_CM_batches"]),
)
tracks_df = pd.read_csv(path + ".tracks.bed", sep="\t", header=None)
content_df = pd.read_csv(path + ".content.txt", sep="\t", header=None)

output_path = os.path.join(
    path_to_output_directory,
    "_".join([dataset_name, "Clomics_CM"]),
)
tracks_df.sort_values([0, 1], ignore_index=True, inplace=True)
new_cm_id = "cm" + tracks_df.index.astype(str)
old_to_new_cm_map = dict(zip(tracks_df.loc[:, 3].to_list(), new_cm_id.to_list()))
tracks_df.sort_values([0, 1], ignore_index=True, inplace=True)
tracks_df.loc[:, 3] = tracks_df.loc[:, 3].map(old_to_new_cm_map)

content_df.loc[:, 0] = content_df.loc[:, 0].map(old_to_new_cm_map)
content_df.loc[:, "index"] = content_df.loc[:, 0].copy()
content_df.set_index("index", inplace=True)
content_df.loc[tracks_df.loc[:, 3], :]

tracks_df.to_csv(output_path + ".tracks.bed", sep="\t", index=False, header=False)
content_df.loc[tracks_df.loc[:, 3], :].to_csv(
    output_path + ".content.txt", sep="\t", index=False, header=False
)
