import argparse
import logging
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import tqdm
import sys

parser = argparse.ArgumentParser(description="Convert CRDs to CMs ©Olga Pushkarev")
parser.add_argument("-d", "--dataset", type=str, help="Input dataset", required=True)
parser.add_argument(
    "-i",
    "--path_to_input_file",
    type=str,
    help="Path to folder with correlation edgelists",
    required=True,
)
parser.add_argument(
    "-o",
    "--output_path",
    type=str,
    help="Output directory for CM tracks and content files",
    required=True,
)
parser.add_argument(
    "-p",
    "--n_processes",
    help="Number of processes to initialize",
    default=1,
    required=False,
)

args = parser.parse_args()

# Arguments
dataset = args.dataset
path_to_input_file = args.path_to_input_file
path_to_output_directory = args.output_path
n_cores = int(args.n_processes)
print("Using ", str(n_cores), " cores")
path_to_output_directory = path_to_output_directory.replace("\\", "/")

logging.basicConfig(
    filename=os.path.join(path_to_output_directory, dataset + "_logfile.txt"),
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    filemode="a",
)
logging.info("########################################")
logging.info("# Convert CRDs to CMs ©Olga Pushkarev #")
logging.info("########################################")
logging.info("\nList of arguments:\n")
logging.info("Input:")
logging.info("\t 1. Input dataset: " + dataset)
logging.info("\t 2. Input path to Clomics .bed file: " + path_to_input_file)
logging.info("\n")
logging.info(
    "\t 1. Output directory for CRDs in CM-like format, i.e. tracks and content files: "
    + path_to_output_directory
)
logging.info("\n")


def split_nested_array(arr, k):
    n = len(arr)
    chunk_size = n // k
    remainder = n % k

    chunks = [
        arr[
            i * chunk_size
            + min(i, remainder) : (i + 1) * chunk_size
            + min(i + 1, remainder)
        ]
        for i in range(k)
    ]
    return chunks


def get_tracks_info(crd_tree_df, crd_nodes_chunk, n_batch):
    with open(
        os.path.join(
            path_to_output_directory, dataset + "_Clomics_CM_batches.tracks.bed"
        ),
        "a+",
    ) as tracks_file, open(
        os.path.join(
            path_to_output_directory, dataset + "_Clomics_CM_batches.content.txt"
        ),
        "a+",
    ) as content_file:

        for index, crd_values in enumerate(crd_nodes_chunk):
            chromosome, start, end, uid, _, parent1, parent2, _ = crd_values
            list_with_marks = find_leafs(
                [parent1, parent2]
            )  # Find marks (leafs of the tree) that belong to one CRD
            # Fill the fields of the .content file
            content_file.write(
                "\t".join(
                    [
                        "cm" + str(index) + "_" + str(n_batch),
                        str(int(len(list_with_marks))),
                        ",".join(list_with_marks),
                    ]
                )
                + "\n"
            )

            # For marks in CM prepare list of peak coordinates (start, end)
            start_end = sorted(
                [
                    (
                        int(
                            crd_tree_df.loc[crd_tree_df["UID"] == mark, "start"].values[
                                0
                            ]
                        ),
                        int(
                            crd_tree_df.loc[crd_tree_df["UID"] == mark, "end"].values[0]
                        ),
                    )
                    for mark in list_with_marks
                ],
                key=lambda x: x[0],
            )
            start = np.array([pair[0] for pair in start_end])
            peak_starts = list(start - min(start_end, key=lambda x: x[0])[0])
            crd_start = min(start_end, key=lambda x: x[0])[0]
            crd_end = max(start_end, key=lambda x: x[1])[1]
            tracks_file.write(
                "\t".join(
                    [
                        str(chromosome),
                        str(crd_start),
                        str(crd_end),
                        "cm" + str(index) + "_" + str(n_batch),
                        str(int(1000)),
                        "+",
                        str(crd_start),
                        str(crd_end),
                        "0,0,0",
                        str(int(len(list_with_marks))),
                        ",".join(map(str, [j - i for i, j in start_end])),
                        ",".join(map(str, peak_starts)),
                        str(
                            1
                            if len(merge_intersecting_intervals(start_end)) == 1
                            else 0
                        ),
                    ]
                )
                + "\n"
            )


# Functions
def find_children(lst):
    """Function that finds children of each node in the input list.
    If the node is not internal (is not epigenetic mark, i.e. H3K4me3, H3K27ac, etc.)
    the functions does not search for its kids and continues to the other nodes.
    If no internal nodes left, the function outputs the list of all children of the nodes in the input list
    """
    children_list = []
    for node in lst:
        if node.startswith("chr"):
            row = crd_tree_df[crd_tree_df["feature_id"] == str(node)]
            if (str(row["first_child"].iloc[0]) == "nan") or (
                str(row["second_child"].iloc[0]) == "nan"
            ):
                children_list.append(str(row["UID"].iloc[0]))
            else:
                children_list.append(str(row["first_child"].iloc[0]))
                children_list.append(str(row["second_child"].iloc[0]))
        else:
            children_list.append(node)
    return children_list


def find_leafs(list_of_parents):
    """Recursive function to find individual peaks of histone marks (a leaf node) that belong to one CRD"""
    if not any(child.startswith("chr") for child in list_of_parents):
        return list_of_parents
    else:
        children = find_children(list_of_parents)
        return find_leafs(children)


def merge_intersecting_intervals(list_of_pairs):
    """For a given list of pairs, the function finds and merges intersecting intervals.
    The output is the list of intervals of length <= length(input list)"""
    merged_intervals = []
    for pair in sorted(list_of_pairs):
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


# Import data from Clomics output
crd_tree_df = pd.read_csv(path_to_input_file, sep="\t", header=None)
crd_tree_df.columns = [
    "chr",
    "start",
    "end",
    "UID",
    "feature_id",
    "first_child",
    "second_child",
    "bool_CRD",
]
crd_tree_df = crd_tree_df.astype({"first_child": str, "second_child": str})

peak_id_to_check = crd_tree_df[crd_tree_df.loc[:, "first_child"] == "nan"].loc[0, "UID"]

if peak_id_to_check.startswith("chr"):
    sys.exit(
        "Error! Rename the peaks to avoid confusion between the names of inner nodes of a tree and its leafs."
    )
else:
    # Select nodes that were marked as CRDs
    crd_nodes = crd_tree_df[crd_tree_df.loc[:, "bool_CRD"] == 1].to_numpy()

    with mp.Pool(processes=int(n_cores)) as pool:
        # starts the sub-processes without blocking
        proc_results = [
            pool.apply_async(
                get_tracks_info, args=(crd_tree_df, crd_nodes_chunk, n_batch)
            )
            for n_batch, crd_nodes_chunk in enumerate(
                split_nested_array(crd_nodes, n_cores)
            )
        ]
        # blocks until all results are fetched
        result_chunks = [r.get() for r in proc_results]
