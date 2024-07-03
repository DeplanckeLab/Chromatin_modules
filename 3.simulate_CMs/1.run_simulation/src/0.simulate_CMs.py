import numpy as np
import os
import pandas as pd
import sys
import warnings
import itertools
from annoy import AnnoyIndex
from multiprocess import Pool
from sklearn import preprocessing

from cm_simulation import Simulation

warnings.filterwarnings("ignore")


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


def simulate_cms(
    method,
    chromosome,
    n_sim_per_cm,
    cm_param_df,
    cm_coordinates_dict,
    ncm_peak_dict,
    k27ac_peak_coord_dict,
    k4me1_peak_coord_dict,
    count_mtx_dict,
    single_peak_list,
):
    simulation_chr = Simulation(
        method=method,
        n_sim_per_cm=n_sim_per_cm,
        chromosome=chromosome,
        cm_df=cm_param_df.loc[cm_param_df.loc[:, "chr"] == chromosome, :],
        cm_peak_coord_dict=cm_coordinates_dict.get(chromosome),
        ncm_peak_coord_dict=ncm_peak_dict.get(chromosome),
        ncm_k27ac_peak_coord_dict=k27ac_peak_coord_dict.get(chromosome),
        ncm_k4me1_peak_coord_dict=k4me1_peak_coord_dict.get(chromosome),
        chr_count_mtx=count_mtx_dict.get(chromosome),
        single_lead_peaks=single_peak_list,
    )
    ref_feature, simulated_cms, simulated_cms_pids, sim_feature_vec = (
        simulation_chr.run_one_simulation()
    )
    return ref_feature, simulated_cms, simulated_cms_pids, sim_feature_vec


def process_chromosome(input_args):
    (
        method,
        chromosome,
        n_sim_per_cm,
        cm_param_df,
        cm_coordinates_dict,
        ncm_peak_dict,
        k27ac_peak_coord_dict,
        k4me1_peak_coord_dict,
        count_mtx_dict,
        single_peak_list,
    ) = input_args
    ref_feature, simulated_cms, simulated_cms_pids, sim_feature_vec = simulate_cms(
        method,
        chromosome,
        n_sim_per_cm,
        cm_param_df,
        cm_coordinates_dict,
        ncm_peak_dict,
        k27ac_peak_coord_dict,
        k4me1_peak_coord_dict,
        count_mtx_dict,
        single_peak_list,
    )
    ref_feature_df = pd.DataFrame.from_dict(ref_feature, orient="index")
    sim_feature_df = pd.concat(
        [
            pd.DataFrame.from_dict(sim_feature, orient="index")
            for _, sim_feature in sim_feature_vec.items()
        ],
        axis=0,
    )
    return chromosome, {
        "ref_feature_df": ref_feature_df,
        "sim_feature_df": sim_feature_df,
        "simulated_cms": simulated_cms,
        "simulated_cms_pids": simulated_cms_pids,
    }


def get_ref_sim_cm_pairs(output_dict_chr, n_closest_cms):
    crd_ref_sim_pairs = []
    used_simulated_cms = set()
    reference_cm_df = output_dict_chr["ref_feature_df"]
    simulated_cm_df = output_dict_chr["sim_feature_df"]
    for crd_id in output_dict_chr["simulated_cms"].keys():
        simulated_cm_df = simulated_cm_df.loc[
            ~simulated_cm_df.index.isin(used_simulated_cms), :
        ]
        reference_cm = reference_cm_df.loc[crd_id, :]

        all_cms = simulated_cm_df.append(reference_cm)
        if all_cms.shape[0] == 1:
            continue

        all_cms = all_cms[
            all_cms.loc[:, "cm_size"] == reference_cm_df.loc[crd_id, "cm_size"]
        ]

        del all_cms["cm_size"]

        min_max_scaler = preprocessing.MinMaxScaler()
        all_cms_np = min_max_scaler.fit_transform(all_cms)
        all_cms_df = pd.DataFrame(
            all_cms_np, columns=all_cms.columns, index=all_cms.index
        )
        f = all_cms_df.shape[1]
        annoy_idx = AnnoyIndex(f, "euclidean")
        annoy_idx.set_seed(42)
        crd_idx_dict = {}
        for i, (s_crd_id, feature_vec) in enumerate(
            all_cms_df.T.to_dict(orient="list").items()
        ):
            annoy_idx.add_item(i, feature_vec)
            crd_idx_dict[i] = s_crd_id
        annoy_idx.build(200)

        closest_sim_cms = [
            crd_idx_dict[el]
            for el in annoy_idx.get_nns_by_item(
                all_cms_df.shape[0] - 1, n_closest_cms + 1
            )
            if crd_idx_dict[el] != crd_id
        ]
        for closest_sim_cm in closest_sim_cms:
            used_simulated_cms.update([closest_sim_cm])
        crd_ref_sim_pairs.append([crd_id, closest_sim_cms])
    return crd_ref_sim_pairs


# --------------

# Arguments
dataset = "MISSING"
path_to_input_directory = "MISSING"
path_to_output_directory = "MISSING"
cm_calling_method = "MISSING"
n_closest_cms = "MISSING"
n_sim_per_cm = "MISSING"
parralelize = "MISSING"
n_cores = "MISSING"

dataset = sys.argv[1]
path_to_input_directory = sys.argv[2]
path_to_output_directory = sys.argv[3]
cm_calling_method = sys.argv[4]
n_closest_cms = int(sys.argv[5])
n_sim_per_cm = int(sys.argv[6])
parralelize = int(sys.argv[7])
n_cores = int(sys.argv[8])
print("#" * 56)
print("# CM simulation ©Olga Pushkarev - List of arguments: #")
print("#" * 56)
print("1. [Input] Input dataset:", dataset)
print("2. [Input] Input path:", path_to_input_directory)
print("3. [Input] cm calling method:", cm_calling_method)
print("4. [Input] Number of closest cms to choose:", n_closest_cms)
print("5. [Input] Number of cms to simulate per reference cm:", n_sim_per_cm)

if parralelize:
    print("-" * 10)
    print("Using ", n_cores, " cores...")
    print("-" * 10)

print("6. [Output] Output path:", path_to_output_directory)
print("#" * 56)


if not os.path.exists(path_to_input_directory):
    sys.exit("Input directory does not exist. Stopping...")
if not os.path.exists(path_to_output_directory):
    os.makedirs(path_to_output_directory)

# Check arguments
if len(sys.argv) != 9:
    sys.exit("Some arguments are MISSING. Stopping...")

# Load the data for simulation:
# Dictionary with count matrices per chromosome
count_mtx_dict = np.load(
    os.path.join(path_to_input_directory, "peak_count_matrix_dict_by_chr.npy"),
    allow_pickle=True,
).item()
# Dataframe with cm features, including number of peaks in a cm, cm length, fraction of H3K27ac, etc
cm_param_df = pd.read_csv(
    os.path.join(path_to_input_directory, "cm_parameter_df.bed"), sep="\t"
)
if not "chr" in str(cm_param_df.iloc[0, 0]):
    cm_param_df.iloc[:, 0] = "chr" + cm_param_df.iloc[:, 0].astype(str)
# Dictionary with cms and the respective peak coordinates
cm_coordinates_dict = np.load(
    os.path.join(path_to_input_directory, "cm_peak_coordinates_dict_by_chr.npy"),
    allow_pickle=True,
).item()
# Dictionary with peak coordinates that were not called as cm-peaks
ncm_peak_dict = np.load(
    os.path.join(path_to_input_directory, "not_cm_peak_coordinates_dict_by_chr.npy"),
    allow_pickle=True,
).item()
if cm_calling_method in ["PHM", "phm"]:
    single_peak_df = pd.read_csv(
        os.path.join(path_to_input_directory, "single_peak_df_for_simulation.txt"),
        sep="\t",
        header=None,
        names=["chromosome", "peak_start", "peak_end", "mark", "peak_id"],
    )
    single_peak_list = single_peak_df.loc[:, "peak_id"].to_list()
else:
    single_peak_list = None

# Split the dictionary with peak coordinates that were not called as cm-peaks
# into histone modification-specific dictionaries
k27ac_peak_coord_dict = {}
k4me1_peak_coord_dict = {}
for chr_id, chr_peak_dict in ncm_peak_dict.items():
    k27ac_peak_dict = {}
    k4me1_peak_dict = {}
    for mark_peak, peak_coord in chr_peak_dict.items():
        if "H3K27ac" in mark_peak:
            k27ac_peak_dict[mark_peak] = peak_coord
        else:
            k4me1_peak_dict[mark_peak] = peak_coord
    k27ac_peak_coord_dict[chr_id] = k27ac_peak_dict
    k4me1_peak_coord_dict[chr_id] = k4me1_peak_dict

chromosomes = sorted(
    set(cm_coordinates_dict.keys()).intersection(set(cm_param_df["chr"]))
)
if parralelize:
    number_of_cores = min(n_cores, len(chromosomes))
    with Pool(int(number_of_cores)) as p:
        pool_outputs = list(
            p.imap(
                process_chromosome,
                zip(
                    itertools.repeat(cm_calling_method),
                    chromosomes,
                    itertools.repeat(n_sim_per_cm),
                    itertools.repeat(cm_param_df),
                    itertools.repeat(cm_coordinates_dict),
                    itertools.repeat(ncm_peak_dict),
                    itertools.repeat(k27ac_peak_coord_dict),
                    itertools.repeat(k4me1_peak_coord_dict),
                    itertools.repeat(count_mtx_dict),
                    itertools.repeat(single_peak_list),
                ),
            )
        )
    ref_sim_cms_by_chr = dict(pool_outputs)
else:
    ref_sim_cms_by_chr = {
        chromosome: process_chromosome(
            (
                cm_calling_method,
                chromosome,
                n_sim_per_cm,
                cm_param_df,
                cm_coordinates_dict,
                ncm_peak_dict,
                k27ac_peak_coord_dict,
                k4me1_peak_coord_dict,
                count_mtx_dict,
                single_peak_list,
            )
        )
        for chromosome in chromosomes
    }

ref_sim_pairs_by_chr = {
    chromosome: get_ref_sim_cm_pairs(chr_df_dict, n_closest_cms)
    for chromosome, chr_df_dict in ref_sim_cms_by_chr.items()
}

np.save(
    os.path.join(
        path_to_output_directory,
        "_".join(
            [
                dataset,
                "ref",
                str(n_closest_cms),
                "sim",
                cm_calling_method,
                "by_chromosome_no_repeating_peaks.npy",
            ]
        ),
    ),
    ref_sim_cms_by_chr,
    allow_pickle=True,
)
np.save(
    os.path.join(
        path_to_output_directory,
        "_".join(
            [
                dataset,
                "ref",
                str(n_closest_cms),
                "sim",
                cm_calling_method,
                "pairs_by_chromosome_no_repeating_peaks.npy",
            ]
        ),
    ),
    ref_sim_pairs_by_chr,
    allow_pickle=True,
)

n_cms = sum(
    [
        len(sim_cm_list[1])
        for _, pairs_dict in ref_sim_pairs_by_chr.items()
        for sim_cm_list in pairs_dict
    ]
)

# Save .content and .tracks files for simulated cms
bed_df = pd.DataFrame(
    columns=[
        "#chr",
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
        "totem_cm",
        "reference_cm",
        "simulated_cm",
    ],
    index=np.arange(n_cms),
)

crd_content = pd.DataFrame(
    columns=["cm_id", "cm_size", "peaks"], index=np.arange(n_cms)
)
counter = 0
for chromosome, ref_sim_pairs in ref_sim_pairs_by_chr.items():
    sim_cm_coordinates = ref_sim_cms_by_chr[chromosome]["simulated_cms"]
    sim_cm_mark_pids = ref_sim_cms_by_chr[chromosome]["simulated_cms_pids"]
    for reference_cm_id, simulated_cm_ids in ref_sim_pairs:
        for simulated_cm_id in simulated_cm_ids:
            sim_cm_peaks = sorted(
                sim_cm_coordinates[simulated_cm_id.split("_")[0]][simulated_cm_id]
            )

            # Fill the fields of the .content file
            crd_content["cm_id"].loc[counter] = "sim_crd" + str(counter)
            crd_content["cm_size"].loc[counter] = int(len(sim_cm_peaks))

            cm_peak_ids = sim_cm_mark_pids[simulated_cm_id.split("_")[0]][
                simulated_cm_id
            ]
            crd_content["peaks"].loc[counter] = ",".join(cm_peak_ids)

            # Fill the fields of the .bed file
            bed_df["#chr"].loc[counter] = str(chromosome)
            bed_df["start"].loc[counter] = min(sim_cm_peaks, key=lambda x: x[0])[0]
            bed_df["end"].loc[counter] = max(sim_cm_peaks, key=lambda x: x[1])[1]
            bed_df["cm_id"].loc[counter] = "sim_crd" + str(counter)
            bed_df["number"].loc[counter] = int(1000)
            bed_df["strain"].loc[counter] = "+"
            bed_df["start_duplicate"].loc[counter] = min(
                sim_cm_peaks, key=lambda x: x[0]
            )[0]
            bed_df["end_duplicate"].loc[counter] = max(
                sim_cm_peaks, key=lambda x: x[1]
            )[1]
            bed_df["numbers"].loc[counter] = "0,0,0"
            bed_df["cm_size"].loc[counter] = int(len(sim_cm_peaks))
            bed_df["peak_length"].loc[counter] = ",".join(
                map(str, [j - i for i, j in sim_cm_peaks])
            )
            start = np.array([pair[0] for pair in sim_cm_peaks])
            peak_starts = list(start - min(sim_cm_peaks, key=lambda x: x[0])[0])
            bed_df["peak_starts"].loc[counter] = ",".join(map(str, peak_starts))
            bed_df["totem_cm"].loc[counter] = (
                1 if len(merge_intersecting_intervals(sim_cm_peaks)) == 1 else 0
            )
            bed_df["reference_cm"].loc[counter] = reference_cm_id
            bed_df["simulated_cm"].loc[counter] = simulated_cm_id
            counter += 1

# Save bed and content files
bed_df.to_csv(
    os.path.join(
        path_to_output_directory,
        "_".join(
            [
                "simulated",
                str(n_closest_cms),
                "closest",
                cm_calling_method,
                dataset,
                "no_repeating_peaks.tracks.bed",
            ]
        ),
    ),
    sep="\t",
    index=False,
    header=False,
)
crd_content.to_csv(
    os.path.join(
        path_to_output_directory,
        "_".join(
            [
                "simulated",
                str(n_closest_cms),
                "closest",
                cm_calling_method,
                dataset,
                "no_repeating_peaks.content.txt",
            ]
        ),
    ),
    sep="\t",
    header=False,
    index=False,
)
