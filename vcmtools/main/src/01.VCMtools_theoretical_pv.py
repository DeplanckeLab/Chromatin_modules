import argparse
import logging
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import utils

parser = argparse.ArgumentParser(
    description="Get peak correlations and theoretical p-values ©Olga Pushkarev"
)
parser.add_argument("-d", "--dataset", type=str, help="Input dataset", required=True)
parser.add_argument(
    "-m",
    "--max_peak_dist",
    help="Maximum distance between peak centers",
    default=0.5,
    required=False,
)
parser.add_argument(
    "-i",
    "--path_to_input_directory",
    type=str,
    help="Path to folder with count matrices",
    default=0,
    required=False,
)
parser.add_argument(
    "-f",
    "--input_files",
    type=str,
    help="Comma separated paths to files with count matrices and mark specification, e.g., H3K4me1:/path/to/count_matrix_1,H3K27ac:/path/to/count_matrix_2,",
    default=0,
    required=False,
)
parser.add_argument(
    "-c", "--chromosomes", type=str, help="Chromosomes to process", required=True
)
parser.add_argument(
    "-n", "--n_cores", help="Number of cores to use", type=int, required=True
)
parser.add_argument(
    "-o",
    "--output_path",
    type=str,
    help="Output directory for correlation values and theoretical p-values",
    required=True,
)

args = parser.parse_args()

# Arguments
dataset = args.dataset
max_peak_dist = float(args.max_peak_dist) * (10**6)
path_to_input_directory = args.path_to_input_directory
input_files = args.input_files
chromosomes_str = args.chromosomes
n_cores = args.n_cores
path_to_output_directory = args.output_path

if "," in chromosomes_str:
    chromosomes = chromosomes_str.split(",")
elif "-" in chromosomes_str:
    chromosomes_lst = [el.replace("chr", "") for el in chromosomes_str.split("-")]
    chromosomes = [
        str(chr_int)
        for chr_int in np.arange(int(chromosomes_lst[0]), int(chromosomes_lst[1]) + 1)
    ]
else:
    chromosomes = [chromosomes_str]

path_to_output_directory = path_to_output_directory.replace("\\", "/")

if not os.path.exists(path_to_output_directory):
    os.makedirs(path_to_output_directory)

logging.basicConfig(
    filename=os.path.join(path_to_output_directory, dataset + "_logfile.txt"),
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    filemode="a",
)
logging.info("##################################################################")
logging.info("# Get peak correlations and theoretical p-values ©Olga Pushkarev #")
logging.info("##################################################################")
logging.info("\nList of arguments:\n")
logging.info("Input:")
logging.info("\t 1. Input dataset: " + dataset)
if path_to_input_directory:
    path_to_input_directory = path_to_input_directory.replace("\\", "/")
    logging.info("\t 2. Path to correlation edgelist(s): " + path_to_input_directory)
elif input_files:
    logging.info("\t 2. Path to input file(s): " + input_files)
logging.info("\t 3. Number of cores to use:" + str(n_cores))
logging.info("\n")
logging.info(
    "\t 1. Output directory for correlation files: " + path_to_output_directory
)
logging.info("\n")

# chromosomes = [str(i) for i in range(1, 23)]
dict_with_count_matrices = {}
if path_to_input_directory:
    if dataset == "iPSCs":
        dict_with_count_matrices[
            "ATAC"
        ] = "ATAC_iPSC_iPSCmeta_only_PCs16_RPKMnorm_regressed_qqnorm_forCRD.bed.gz"
    elif dataset == "Gaffney":
        dict_with_count_matrices[
            "ATAC"
        ] = "ATAC_Gaffney_gafmeta_only_16PCs_RPKMnorm_regressed_qqnorm_forCRD.bed.gz"
    else:
        for file_name in os.listdir(os.path.join(path_to_input_directory, dataset)):
            if file_name.endswith("_forCRD.bed.gz"):
                dict_with_count_matrices[file_name.split("_")[0]] = file_name
elif input_files:
    input_files_list = input_files.split(",")
    for input_file in input_files_list:
        dict_with_count_matrices[input_file.split(":")[0]] = input_file.split(":")[1]

dict_with_dataset_dfs = {}
for mark, mark_count_mtx_path in dict_with_count_matrices.items():
    if path_to_input_directory:
        if (dataset == "iPSCs") or (dataset == "Gaffney"):
            df = pd.read_csv(
                os.path.join(path_to_input_directory, mark_count_mtx_path), sep="\t"
            )
            df["#Chr"] = df["#Chr"].str.replace("chr", "")
        else:
            df = pd.read_csv(
                os.path.join(path_to_input_directory, dataset, mark_count_mtx_path),
                sep="\t",
            )
    elif input_files:
        df = pd.read_csv(mark_count_mtx_path, sep="\t")
    df.rename({"#chr": "#Chr", "id": "pid", "dummy": "gid"}, axis=1, inplace=True)
    if not "chr" in df.iloc[0, 3]:
        df["pid"] = "chr" + df["pid"]
    if not mark in df.iloc[0, 3]:
        df["pid"] = mark + ":" + df["pid"]
    df = df.set_index("pid")
    df["#Chr"] = df["#Chr"].astype(str)
    df.drop(["start", "end", "gid", "strand"], axis=1, inplace=True)
    dict_with_dataset_dfs[mark] = df

counts_for_chr_dict = {}
for chromosome in chromosomes:
    all_data_df_chr = pd.concat(
        [
            mark_df[mark_df.loc[:, "#Chr"] == chromosome]
            for mark, mark_df in dict_with_dataset_dfs.items()
        ],
        axis=0,
        join="outer",
        sort=False,
    )
    if "#Chr" in all_data_df_chr.columns:
        all_data_df_chr.drop(["#Chr"], axis=1, inplace=True)
    counts_for_chr_dict[chromosome] = all_data_df_chr.T.to_dict("list")

with mp.Pool(processes=int(n_cores)) as pool:
    # starts the sub-processes without blocking
    proc_results = [
        pool.apply_async(
            utils.process_different_marks_pv,
            args=(path_to_output_directory, dataset, chunk, max_peak_dist),
        )
        for chromosome, chunk in counts_for_chr_dict.items()
    ]
    # blocks until all results are fetched
    result_chunks = [r.get() for r in proc_results]
