import argparse
import os
import pandas as pd
import sys
import warnings
import utils

warnings.filterwarnings("ignore")


def is_unique(s):
    a = s.to_numpy()
    return (a[0] == a).all()


parser = argparse.ArgumentParser(
    description="Match simulated CMs with mapped (reference) ©Olga Pushkarev"
)
# Output path
parser.add_argument(
    "-i",
    "--input_path",
    type=str,
    help="Input directory with simulated CMs tracks and content files",
    required=True,
)
parser.add_argument("-d", "--dataset", type=str, help="Input dataset", required=True)
parser.add_argument(
    "-m",
    "--method",
    type=str,
    help="CM mapping method",
    required=True,
)
# Path to mapped CMs
parser.add_argument(
    "-cmp",
    "--core_cm_path",
    help="Path to mapped CMs",
    type=str,
    required=True,
)

# CM simulation parameter
parser.add_argument(
    "-n",
    "--n_closest_cms",
    help="Number of closest CMs to choose; parameter for simulated CMs",
    required=True,
)
parser.add_argument(
    "-k",
    "--k_sim_per_cm",
    help="Number of CMs to simulate per reference CMs",
    required=True,
)
# PHM parameters
parser.add_argument(
    "-pp",
    "--pp_threshold",
    help="Posterior probability of causal intercations for PHM",
    default="0.8",
    required=False,
)
# VCMtools parameters
parser.add_argument(
    "-vw",
    "--vcm_window",
    help="cis- window used in VCMtools",
    default="0.5Mb",
    type=str,
    required=False,
)
parser.add_argument(
    "-pv",
    "--pv_threshold",
    help="p-value threshold used in VCMtools",
    default=0.001,
    required=False,
)
# Clomics parameters
parser.add_argument(
    "-np",
    "--n_peaks",
    help="Number of peaks used in Clomics",
    default="n_peaks_200",
    type=str,
    required=False,
)
parser.add_argument(
    "-bg",
    "--bg_threshold",
    help="Background correlation threshold used in Clomics",
    default="bg_threshold_3",
    type=str,
    required=False,
)
# # Output path
# parser.add_argument(
#     "-o",
#     "--output_path",
#     type=str,
#     help="Output directory for formated CM tracks and content files",
#     required=True,
# )

args = parser.parse_args()
path_to_input_directory = args.input_path

core_cm_path = args.core_cm_path
dataset = args.dataset
method = args.method
n_closest_cms = str(args.n_closest_cms)
k_sim_per_cm = str(args.k_sim_per_cm)
pp_threshold = str(args.pp_threshold)
vcm_window = str(args.vcm_window)
pv_threshold = str(args.pv_threshold)
n_peaks = str(args.n_peaks)
bg_threshold = str(args.bg_threshold)

# path_to_output_directory = args.output_path


path_to_input_directory = path_to_input_directory.replace("\\", "/")
if not os.path.exists(path_to_input_directory):
    sys.exit("Input directory does not exist. Stopping...")

path_to_output_directory = os.path.join(
    path_to_input_directory, method, "all_simulated_cms"
)
if not os.path.exists(path_to_output_directory):
    os.makedirs(path_to_output_directory)

tracks_df_bed_ref, content_df_bed_ref = utils.get_cm_tracks_content(
    core_cm_path,
    method,
    dataset,
    pp_threshold=pp_threshold,
    vcm_window=vcm_window,
    pv_threshold=pv_threshold,
    n_peaks=n_peaks,
    bg_threshold=bg_threshold,
)
tracks_df_bed_ref.loc[:, 3] = tracks_df_bed_ref.loc[:, 3].str.replace("_", "~")
content_df_bed_ref.loc[:, 0] = content_df_bed_ref.loc[:, 0].str.replace("_", "~")

tracks_df_bed_sim, content_df_bed_sim = utils.get_simulated_cm_tracks_content(
    path_to_input_directory,
    dataset,
    method,
    n_closest_cms=n_closest_cms,
    k_sim_per_cm=k_sim_per_cm,
)

cm_overlap = set(tracks_df_bed_ref.loc[:, 3]).intersection(
    set(tracks_df_bed_sim.loc[:, 13])
)
tracks_ref = tracks_df_bed_ref.loc[tracks_df_bed_ref.loc[:, 3].isin(cm_overlap), :]
content_ref = content_df_bed_ref.loc[content_df_bed_ref.loc[:, 0].isin(cm_overlap), :]
tracks_sim = tracks_df_bed_sim.loc[tracks_df_bed_sim.loc[:, 13].isin(cm_overlap), :]
content_sim = content_df_bed_sim.loc[
    content_df_bed_sim.loc[:, 0].isin(tracks_sim.loc[:, 3]), :
]

if tracks_ref.shape[0] == content_ref.shape[0]:
    tracks_ref.to_csv(
        os.path.join(
            path_to_output_directory,
            "_".join([dataset, "ref", method, "no_repeating_peaks_matched.tracks.bed"]),
        ),
        sep="\t",
        index=False,
        header=False,
    )

    content_ref.to_csv(
        os.path.join(
            path_to_output_directory,
            "_".join(
                [dataset, "ref", method, "no_repeating_peaks_matched.content.txt"]
            ),
        ),
        sep="\t",
        index=False,
        header=False,
    )

if tracks_sim.shape[0] == content_sim.shape[0]:
    tracks_sim.to_csv(
        os.path.join(
            path_to_output_directory,
            "_".join([dataset, "sim", method, "no_repeating_peaks_matched.tracks.bed"]),
        ),
        sep="\t",
        index=False,
        header=False,
    )

    content_sim.to_csv(
        os.path.join(
            path_to_output_directory,
            "_".join(
                [dataset, "sim", method, "no_repeating_peaks_matched.content.txt"]
            ),
        ),
        sep="\t",
        index=False,
        header=False,
    )
else:
    print("Shapes are not okay!")

# # Save filtered tracks and content
# if not os.path.exists(os.path.join(path_to_output_directory, "all_simulated_cms")):
#     os.makedirs(os.path.join(path_to_output_directory, "all_simulated_cms"))

# tracks_df.to_csv(
#     os.path.join(
#         path_to_output_directory,
#         "_".join([dataset, "sim", method, "no_repeating_peaks.tracks.bed"]),
#     ),
#     sep="\t",
#     header=False,
#     index=False,
# )
# cm_content.to_csv(
#     os.path.join(
#         path_to_output_directory,
#         "_".join([dataset, "sim", method, "no_repeating_peaks.content.txt"]),
#     ),
#     sep="\t",
#     header=False,
#     index=False,
# )
