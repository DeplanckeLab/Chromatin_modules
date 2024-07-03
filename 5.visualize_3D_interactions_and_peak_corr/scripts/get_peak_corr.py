import argparse
import logging
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import scipy.stats as stats

parser = argparse.ArgumentParser(
    description="Get peak correlations and theoretical p-values ©Olga Pushkarev"
)
parser.add_argument("-d", "--dataset", type=str, help="Input dataset", required=True)
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
logging.info("\t 2. Path to input file(s): " + input_files)
logging.info("\t 3. Number of cores to use:" + str(n_cores))
logging.info("\n")
logging.info(
    "\t 1. Output directory for correlation files: " + path_to_output_directory
)
logging.info("\n")


def process_different_marks_pv(
    path_to_output_directory, dataset, dict_for_chr, max_dist=5 * (10**6)
):
    with open(
        os.path.join(
            path_to_output_directory,
            dataset + "_theoretical_corr_with_p_values.bed",
        ),
        "a+",
    ) as file_with_corr:
        for peak1, peak1_values in dict_for_chr.items():
            peak1_variables = peak1.split(":")
            mark1 = peak1_variables[0]
            chromosome_1 = peak1_variables[1]
            if not "chr" in str(chromosome_1):
                chromosome_1 = "chr" + str(chromosome_1)
            peak1_start = int(peak1_variables[2])
            peak1_end = int(peak1_variables[3])
            peak1_center = peak1_start + (peak1_end - peak1_start) / 2
            for peak2, peak2_values in dict_for_chr.items():
                peak2_variables = peak2.split(":")
                mark2 = peak2_variables[0]
                chromosome_2 = peak2_variables[1]
                if not "chr" in str(chromosome_2):
                    chromosome_2 = "chr" + str(chromosome_2)
                peak2_start = int(peak2_variables[2])
                peak2_end = int(peak2_variables[3])
                peak2_center = peak2_start + (peak2_end - peak2_start) / 2
                try:
                    if (peak1 != peak2) and (
                        abs(peak1_center - peak2_center) <= max_dist
                    ):
                        x = np.array(peak1_values)
                        y = np.array(peak2_values)
                        mask = ~np.isnan(x) & ~np.isnan(y)
                        x = x[mask]
                        y = y[mask]
                        theoretical_correlation, p_value = stats.pearsonr(x, y)
                        file_with_corr.write(
                            "\t".join(
                                [
                                    chromosome_1,
                                    str(peak1_start),
                                    str(peak1_end),
                                    peak1,
                                    chromosome_2,
                                    str(peak2_start),
                                    str(peak2_end),
                                    peak2,
                                    str(theoretical_correlation),
                                    #                                                         str(p_value)
                                ]
                            )
                            + "\n"
                        )
                except:
                    continue


input_files_list = input_files.split(",")
dict_with_count_matrices = {}
for input_file in input_files_list:
    dict_with_count_matrices[input_file.split(":")[0]] = input_file.split(":")[1]

dict_with_dataset_dfs = {}
for mark, mark_count_mtx_path in dict_with_count_matrices.items():
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
            mark_df[mark_df["#Chr"] == chromosome]
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
            process_different_marks_pv, args=(path_to_output_directory, dataset, chunk)
        )
        for chromosome, chunk in counts_for_chr_dict.items()
    ]
    # blocks until all results are fetched
    result_chunks = [r.get() for r in proc_results]
