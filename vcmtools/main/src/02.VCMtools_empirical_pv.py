import argparse
import logging
import os
import utils

parser = argparse.ArgumentParser(
    description="Build VCMs from correlation edgelist ©Olga Pushkarev"
)
parser.add_argument("-d", "--dataset", type=str, help="Input dataset", required=True)
parser.add_argument(
    "-i",
    "--path_to_input_directory",
    type=str,
    help="Path to folder with theoretical p-values for correlations",
    required=True,
)
parser.add_argument(
    "-o",
    "--output_path",
    type=str,
    help="Output directory for correlation values and empirical p-values",
    required=True,
)

args = parser.parse_args()

# Arguments
dataset = args.dataset
path_to_input_directory = args.path_to_input_directory
path_to_output_directory = args.output_path

path_to_input_directory = path_to_input_directory.replace("\\", "/")
path_to_output_directory = path_to_output_directory.replace("\\", "/")

if not os.path.exists(path_to_output_directory):
    os.makedirs(path_to_output_directory)

logging.basicConfig(
    filename=os.path.join(path_to_output_directory, dataset + "_logfile.txt"),
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    filemode="a",
)
logging.info("###########################################################")
logging.info("# Get empirical p-values for correlations ©Olga Pushkarev #")
logging.info("###########################################################")
logging.info("\nList of arguments:\n")
logging.info("Input:")
logging.info("\t 1. Input dataset: " + dataset)
logging.info("\t 2. Path to correlation edgelist(s): " + path_to_input_directory)
logging.info("\n")
logging.info(
    "\t 1. Output directory for correlation files: " + path_to_output_directory
)
logging.info("\n")

input_path = os.path.join(
    path_to_input_directory,
    dataset + "_all_marks_VCM_theoretical_corr_with_p_values.txt",
)

utils.empirical_pvalue_for_corr(dataset, input_path, path_to_output_directory)
