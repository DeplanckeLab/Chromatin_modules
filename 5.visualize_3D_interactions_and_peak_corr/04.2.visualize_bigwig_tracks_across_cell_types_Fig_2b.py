import pandas as pd
import os
import numpy as np
import random
import sys
import tqdm.auto as tqdm

random.seed(42)
sys.path.append(
    "/data/pushkare/computational_paper/06.peak_based_3D_interactions_and_correlations"
)
import peak_matrix_visualization as pmv

# %load_ext autoreload
# %autoreload 2

# # +
tracks_path = "/data/pushkare/computational_paper/01.mapping_of_chromatin_modules"
cm_peak_path = "/data/pushkare/computational_paper/03.peaks_in_chromatin_modules"

corr_mtx_path = os.path.join(
    "/data/pushkare/computational_paper/06.peak_based_3D_interactions_and_correlations",
    "peak_correlations",
)
core_path = "/data/pushkare/computational_paper/06.peak_based_3D_interactions_and_correlations/TADs_CTCF_CM_tracks"
aCM_scores_path = os.path.join(tracks_path, "aCM_scores")

# # +
datasets = ["FIB", "LCL", "Monocytes", "Neutrophils", "Tcells"]
marks = "all_marks"
pp_threshold = 0.8

color_dict = {
    "LCL": "#8DA0CB",
    "FIB": "#FC8D62",
    "Monocytes": "#66C2A5",
    "Neutrophils": "#E78AC3",
    "Tcells": "#A6D854",
}
# -

# ### Get gene info

# Get genes and exon coordinates per gene
exons_per_chr_dict = np.load(
    "/data/pushkare/genome/hg19/all_gene_exon_coordinates_per_chr.npy",
    allow_pickle=True,
).item()
protein_coding_lincRNA_genes = pd.read_csv(
    os.path.join(
        "/data/pushkare/computational_paper/04.functional_annotation_of_chromatin_modules",
        "hg19_gene_annotations",
        "protein_coding_linc_RNA_gene_coordinates_and_ensembl_hgnc_ids_with_chr.bed",
    ),
    sep="\t",
    header=None,
    names=[
        "chr",
        "start",
        "end",
        "strand",
        "ensembl_gene_id",
        "gene_symbol",
        "gene_type",
        "category",
    ],
)
# Get dictionary with ENSEMBL to gene symbol mapping
ens_id_gene_symbol_mapping = dict(
    zip(
        protein_coding_lincRNA_genes.loc[:, "ensembl_gene_id"],
        protein_coding_lincRNA_genes.loc[:, "gene_symbol"].str.replace(" ", ""),
    )
)
gene_symbol_coord_mapping = dict(
    zip(
        protein_coding_lincRNA_genes.loc[:, "gene_symbol"].str.replace(" ", ""),
        list(
            zip(
                protein_coding_lincRNA_genes.loc[:, "start"],
                protein_coding_lincRNA_genes.loc[:, "end"],
            )
        ),
    )
)
# Extract only protein coding genes
protein_coding_genes = protein_coding_lincRNA_genes.loc[
    protein_coding_lincRNA_genes.loc[:, "gene_type"] == "protein_coding", :
]


target_regions = pd.DataFrame(
    [
        #         ['chr2', 110499915, 111452403, 'chr2_110499915_111452403_LIMS3-RGPD5', '', '', ''],
        #         ['chr7', 102119310, 102347034, 'chr7_102119310_102347034_RASA4', '', '', ''],
        #         ['chr7', 74055331, 75017739, 'chr7_74055331_75017739_GTF2I', '', '', ''],
        #         ['chr15', 43824815, 44090013, 'chr15_43824815_44090013_CATSPER2', '', '', ''],
        #         ['chr15', 20234201, 22689873, 'chr15_20234201_22689873_OR4M2', '', '', ''],
        #         ['chr17', 44349964, 44836662, 'chr17_44349964_44836662_NSF', '', '', ''],
        #         ['chr22', 24296329, 24329551, 'chr22_24296329_24329551_DDTL', '', '', ''],
        #         ['chr8', 145143491, 145525516, 'chr8_145143491_145525516_HGH1', '', '', ''],
        #         ['chr3', 192350398, 192800636, 'chr3_192350398_192800636_MB21D2', '', '', ''],
        #         ['chr4', 108897182, 109147810, 'chr4_108897182_109147810_LEF1', '', '', ''],
        #         ['chr12', 8749853, 8819628, 'chr12_8749853_8819628_AICDA', 'LCL', '', ''],
        #         ['chr11', 33858110, 34008587, 'chr11_33858110_34008587_LMO2', 'LCL', '', ''],
        #         ['chr18', 12750083, 12935149, 'chr18_12750083_12935149_PTPN2', 'LCL', '', ''],
        #         ['chr20', 44709762, 44783632, 'chr20_44709762_44783632_CD40', 'Monocytes', '', ''],
        #         ['chr6', 131887793, 131912924, 'chr6_131887793_131912924_ARG1', 'Neutrophils', '', ''],
        #         ['chr12', 7502452, 7749791, 'chr12_7502452_7749791_CD163', 'Neutrophils', '', ''],
        #         ['chr6', 26400028, 26440047, 'chr6_26400028_26440047_BTN3A1', 'Tcells', '', ''],
        #         ['chr17', 6894279, 6994213, 'chr17_6894279_6994213_BCL6B', 'Tcells', '', ''],
        #         ['chr3', 3099942, 3230733, 'chr3_3099942_3230733_CRBN-IL5RA', 'LCL_Tcells'],
        #         ['chr5', 148736958, 148785027, 'chr5_148736958_148785027_IL17B', 'Tcells', '', ''],
        #         ['chr1', 161538840, 161648886, 'chr1_161538840_161648886_FCGR', '', '', ''],
        #         ['chr1', 25579245, 25785656, 'chr1_25579245_25785656_RHD', '', '', '']
        #         ['chr3', 3106257, 3170695, 'chr3_3106257_3170695_IL5RA', True, '', ''],
        ["chr1", 25579793, 25670256, "chr1_25579793_25670256_RHD-main", True, "", ""],
#         ["chr1", 25770567, 25781098, "chr1_25770567_25781098_RHD-right", False, "", ""],
        ["chr1", 25707947, 25779658, "chr1_25707947_25779658_RHD-right", False, "", ""],
        ["chr2", 61016393, 61231263, "chr2_61016393_61231263_REL", True, "", ""],
        [
            "chr3",
            192505049,
            192645043,
            "chr3_192505049_192645043_MB21D2-main",
            False,
            "",
            "",
        ],
        [
            "chr3",
            192392753,
            192397559,
            "chr3_192392753_192397559_MB21D2-left",
            True,
            "",
            "",
        ],
        [
            "chr3",
            192743619,
            192748425,
            "chr3_192743619_192748425_MB21D2-right",
            False,
            "",
            "",
        ],
        #         ['chr6', 21286083, 23495682, 'chr6_21286083_23495682_SOX4', True, '', '']
        #         ['chr2', 191867259, 192023330, 'chr2_191867259_192023330_STAT4', True, '', ''],
        #         ['chr2', 213740156, 214123694, 'chr2_213740156_214123694_IKZF2', True, '', ''],
        #         ['chr6', 21552226, 22350816, 'chr6_21552226_22350816_SOX4', True, '', ''],
        #         ['chr8', 59613936, 60047497, 'chr6_59613936_60047497_TOX', True, '', ''],
        #         ['chr2', 60636893, 60786250, 'chr2_60636893_60786250_BCL11A', True, '', '']
        #         ['chr12', 101784301, 101950095, 'chr12_101784301_101950095_SPIC', True, '', '']
        #         ['chr14', 96159302, 96197028, 'chr14_96159302_96197028_TCL1A', True, '', '']
    ],
    columns=["chr_id", "start", "end", "pid", "cm_id", "rs_id", "LD_rs_ids"],
)

regions_for_max = {
    "RHD": ["chr1", 25579793, 25779658],
    "MB21D2": ["chr3", 192392753, 192748425],
    "REL": ["chr2", 61016393, 61231263],
    #     'TCL1A': ['chr14', 96159302, 96197028]
    #     'IL5RA': ['chr3', 3106257, 3170695],
    #     'SOX4': ['chr6', 21552226, 22350816],
    #     'STAT4': ['chr2', 191867259, 192023330],
    #     'IKZF2': ['chr2', 213740156, 214123694],
    #     'TOX': ['chr8', 59613936, 60047497],
    #     'BCL11A': ['chr2', 60636893, 60786250]
    #     'SPIC': ['chr12', 101784301, 101950095]
}

path_to_bw_dict = {
    "FIB": "/home/pushkare/NAS1/pushkare/delaneau_bigwig/bigwig",
    "LCL": "/home/pushkare/NAS1/pushkare/delaneau_bigwig/bigwig",
    "Monocytes": "/home/pushkare/NAS1/pushkare/egad2672_monocytes/bigwig",
    "Neutrophils": "/home/pushkare/NAS1/pushkare/egad270_neutrophils/bigwig",
    "Tcells": "/home/pushkare/NAS1/pushkare/egad2673_tcells/bigwig",
}

# region
# Select a sample per call type based on pseudo aCM score
# calculated for peaks of the largest CM in the locus
# in a cell type-specific fashion

max_min_samples_per_locus = {}
paths_to_samples_per_locus = {}
for gene_symbol in ["MB21D2", "RHD", "REL"]:
    max_min_samples_per_ds = {}
    paths_to_samples_per_ds = {}
    for dataset in datasets:
        pseudo_aCM_mtx = pd.read_csv(
            os.path.join(
                "/data/pushkare/computational_paper",
                "06.peak_based_3D_interactions_and_correlations",
                "Fig_2b.pseudo_aCM_mtx_for_largest_CM_in_locus",
                "_".join([gene_symbol, dataset, "pseudo_aCM_matrix.txt"]),
            ),
            sep="\t",
        ).T
        pseudo_aCM_mtx.columns = ["aCM_score"]
        min_sample = pseudo_aCM_mtx.loc[:, "aCM_score"].idxmin()
        max_sample = pseudo_aCM_mtx.loc[:, "aCM_score"].idxmax()
        
        max_min_samples_per_ds[dataset] = {
            "min_sample": min_sample,
            "max_sample": max_sample
        }
        
        path_to_bw = path_to_bw_dict.get(dataset)
        
        paths_to_samples = []
        for sample in [min_sample, max_sample]:
            for f in os.listdir(path_to_bw):
                if 'H3K27ac' in f:
                    condition_k27ac = 1
                    h_mark = 'H3K27ac'
                elif 'H3K4me1' in f:
                    condition_k4me1 = 1
                    h_mark = 'H3K4me1'
                else:
                    continue
                if dataset in ['LCL', 'FIB']:
                    if (sample in f) and (condition_k27ac | condition_k4me1) and (dataset in f):
                        paths_to_samples.append([sample, h_mark, f])
                else:
                    if (sample in f) and (condition_k27ac | condition_k4me1):
                        paths_to_samples.append([sample, h_mark, f])
        paths_to_samples_per_ds[dataset] = paths_to_samples
        
    max_min_samples_per_locus[gene_symbol] = max_min_samples_per_ds
    paths_to_samples_per_locus[gene_symbol] = paths_to_samples_per_ds
# endregion

max_min_samples_per_locus.get('REL').get('LCL')

paths_to_samples_per_locus.get('REL').get('Tcells')

# region

# # Select a sample per cell type
# samples_per_ds = {}
# for dataset, path_to_bw in path_to_bw_dict.items():
#     if dataset in ['LCL', 'FIB']:
#         separator = '_'
#         samples_list = [
#             f.split(separator)[0]
#             for f in os.listdir(path_to_bw)
#             if f.startswith('HG') & (dataset in f)
#         ]
#     else:
#         separator = '.'
#         samples_list = [
#             f.split(separator)[0]
#             for f in os.listdir(path_to_bw)
#         ]
#     samples_per_ds[dataset] = list(set(samples_list))
# endregion

# region
# paths_to_samples_per_ds = {}
# for dataset, path_to_bw in path_to_bw_dict.items():
# #     if dataset == 'LCL':
# #         samples = ['HG00249']
# #     elif dataset == 'FIB':
# #         samples = ['UC1000']
# #     elif dataset == 'Monocytes':
# #         samples = ['S006UKH5']
# #     elif dataset == 'Tcells':
# #         samples = ['S00FULH3'] # ['S00KEX', 'S00KKL', 'S00PDF', 'S0130A', 'S006UK', 'S01342']#['S00FULH3']
# #     elif dataset == 'Neutrophils':
# #         samples = ['S003MBH2']
# #     else:
# #         continue

#     paths_to_samples = []
#     for sample in samples:
#         for f in os.listdir(path_to_bw):
#             if 'H3K27ac' in f:
#                 condition_k27ac = 1
#                 h_mark = 'H3K27ac'
#             elif 'H3K4me1' in f:
#                 condition_k4me1 = 1
#                 h_mark = 'H3K4me1'
#             else:
#                 continue
#             if dataset in ['LCL', 'FIB']:
#                 if (sample in f) and (condition_k27ac | condition_k4me1) and (dataset in f):
#                     paths_to_samples.append([sample, h_mark, f])
#             else:
#                 if (sample in f) and (condition_k27ac | condition_k4me1):
#                     paths_to_samples.append([sample, h_mark, f])
#     paths_to_samples_per_ds[dataset] = paths_to_samples
# endregion

dataset = 'Tcells'
cm_id = "crd1139"
aCM_scores_df = pd.read_csv(
    os.path.join(
        aCM_scores_path,
        marks,
        'Clomics',
        dataset,
        'aCM_matrix',
        'sign_corrected_aCM_matrix.bed'
    ),
    sep='\t',
    index_col=0
)
sorted_aCM_scores = aCM_scores_df.loc[
    cm_id, :
].to_frame().dropna().sort_values(cm_id)
# min_samples = list(sorted_aCM_scores.iloc[:3, :].index)
# max_samples = list(sorted_aCM_scores.iloc[-2:, :].index) + [sorted_aCM_scores.iloc[-3, :].name]

color_ = None
for dataset in tqdm.tqdm(datasets):
    #     if dataset == 'FIB':
    #         continue
    query = zip(
        target_regions["chr_id"],
        target_regions["start"],
        target_regions["end"],
        target_regions["pid"],
        target_regions["cm_id"],
        target_regions["rs_id"],
        target_regions["LD_rs_ids"],
    )
    # sample_info_list = paths_to_rnd_samples_per_ds.get(dataset)
    for chromosome, region_start, region_end, pid, is_left, _, _ in query:
        gene_name = pid.split("_")[-1]
        gene_symbol = gene_name.split("-")[0]
        
        min_max_samples = max_min_samples_per_locus.get(gene_symbol).get(dataset)
        min_sample = min_max_samples.get("min_sample")
        max_sample = min_max_samples.get("max_sample")
        
        bw_files_for_samples = paths_to_samples_per_locus.get(gene_symbol).get(dataset)

        full_region = regions_for_max.get(gene_name.split("-")[0])
        full_region_chr, full_region_start, full_region_end = full_region

        k27ac_maxes = []
        k4me1_maxes = []
        for sample_id, mark, file_name in bw_files_for_samples:
            y_max_ = pmv.get_max_signal_bigwig(
                full_region_chr,
                full_region_start,
                full_region_end,
                bw_path=path_to_bw_dict.get(dataset),
                file_name=file_name,
            )
            if mark == "H3K27ac":
                k27ac_maxes.append(y_max_)
            elif mark == "H3K4me1":
                k4me1_maxes.append(y_max_)

        k27ac_max = np.max(k27ac_maxes)
        k4me1_max = np.max(k4me1_maxes)

        y_max_per_mark = {}
        for sample_id, mark, file_name in bw_files_for_samples:
            if mark == "H3K27ac":
                y_max_per_mark[sample_id + "_" + mark] = k27ac_max
            elif mark == "H3K4me1":
                y_max_per_mark[sample_id + "_" + mark] = k4me1_max

        fig_output_path = os.path.join(core_path, "plots", dataset, gene_name)
        if not os.path.exists(fig_output_path):
            os.makedirs(fig_output_path)
        if "main" in pid:
            if "RHD" in pid:
                fig_width = 7.5
            else:
                fig_width = 5
        elif ("left" in pid) or ("right" in pid):
            fig_width = 2.5
        else:
            fig_width = 7.5  # 10
        
        if "RHD" in pid:
            bin_factor = 0.075
        else:
            bin_factor = 0.075
        if dataset != "CLL":
            for sample_id, mark, file_name in bw_files_for_samples:
                if sample_id == min_sample:
                    if mark == "H3K27ac":
                        color_ = "#98B312"
#                         color_ = "#A1CE00"
                    elif mark == "H3K4me1":
                        color_ = "#D78EFF"
                        
                elif sample_id == max_sample:
                    if mark == "H3K27ac":
                        color_ = "#4E5922"
#                         color_ = "#6D7A3E"
                    elif mark == "H3K4me1":
                        color_ = "#4A4A96"
                else:
                    print('weird stuff')
                    
                pmv.plot_bigwig_profile(
                    dataset,
                    chromosome,
                    region_start,
                    region_end,
                    sample_id,
                    file_name=file_name,
                    mark=mark,
                    bw_path=path_to_bw_dict.get(dataset),
                    bin_factor=bin_factor,
                    fig_width=fig_width,
                    save_fig=True,
                    output_path=fig_output_path,
                    pid=pid,
                    y_max=y_max_per_mark.get(sample_id + "_" + mark),
                    is_left=is_left,
                    color=color_
                )

# ### Plot CV and Standard deviation profiles

# region
# paths_to_samples_per_ds = {}
# for dataset, path_to_bw in path_to_bw_dict.items():
#     samples = samples_per_ds.get(dataset)
#     k27ac_samples_f_names = []
#     k4me1_samples_f_names = []
#     for f in os.listdir(path_to_bw):
#         if "H3K27ac" in f:
#             k27ac_samples_f_names.append(f)
#         elif "H3K4me1" in f:
#             k4me1_samples_f_names.append(f)
#     paths_to_samples_per_ds[dataset] = {
#         "H3K27ac": k27ac_samples_f_names,
#         "H3K4me1": k4me1_samples_f_names,
#     }
# endregion

# region
# query = zip(
#     target_regions["chr_id"],
#     target_regions["start"],
#     target_regions["end"],
#     target_regions["pid"],
#     target_regions["cm_id"],
#     target_regions["rs_id"],
#     target_regions["LD_rs_ids"],
# )
# for cv_std in ["CV", "std"]:
#     for chromosome, region_start, region_end, pid, _, _, _ in query:
#         gene_name = pid.split("_")[-1]
#         fig_output_path = os.path.join(core_path, "plots", "all_datasets")
#         if not os.path.exists(fig_output_path):
#             os.makedirs(fig_output_path)
#         for mark in ["H3K27ac", "H3K4me1"]:
#             fig, ax = plt.subplots(1, 1, figsize=(20, 1))
#             label_patches = []
#             for dataset in datasets:
#                 if dataset == "FIB":
#                     continue
#                 pmv.get_bp_profile_from_bigwig(
#                     chromosome,
#                     region_start,
# region _end,
#                     file_name_list=paths_to_samples_per_ds.get(dataset).get(mark),
#                     #                 mark=mark,
#                     bw_path=path_to_bw_dict.get(dataset),
#                     color=color_dict.get(dataset),
#                     ax=ax,
#                     method=cv_std,
#                 )
#                 label_patches.append(
#                     matplotlib.patches.Patch(
#                         color=color_dict.get(dataset), alpha=1, label=dataset
#                     )
#                 )
#             if mark == "H3K27ac":
#                 ax.legend(
#                     handles=label_patches,
#                     loc="lower center",
#                     bbox_to_anchor=(0.5, -1),
#                     fancybox=False,
#                     shadow=False,
#                     ncol=4,
#                 )
#             plt.suptitle(gene_name + ", " + mark, x=0.6, y=1.05)
#             plt.savefig(
#                 os.path.join(
#                     fig_output_path,
#                     "_".join(
#                         ["all_datasets", pid, mark, "bigwig_bp", cv_std, "profile.pdf"]
#                     ),
#                 ),
#                 bbox_inches="tight",
#                 edgecolor=None,
#                 facecolor=None,
#                 transparent=True,
#                 dpi=300,
#             )
# endregion
