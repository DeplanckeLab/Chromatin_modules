# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.7
#   kernelspec:
#     display_name: Python v3.9
#     language: python
#     name: python3
# ---

import os
import pandas as pd

dataset = "LCL"

for bg_threshold in bg_thresholds:
    path = os.path.join(
        core_path,
        window,
        bg_threshold,
        dataset,
        "_".join([dataset, "Clomics_CM_batches"]),
    )
    tracks_df = pd.read_csv(path + ".tracks.bed", sep="\t", header=None)
    content_df = pd.read_csv(path + ".content.txt", sep="\t", header=None)

    output_path = os.path.join(
        core_path, window, bg_threshold, dataset, "_".join([dataset, "Clomics_CM"])
    )
    tracks_df.sort_values([0, 1], ignore_index=True, inplace=True)
    new_cm_id = "cm" + tracks_df.index.astype(str)
    old_to_new_cm_map = dict(
        zip(tracks_df.loc[:, 3].to_list(), new_cm_id.to_list())
    )
    tracks_df.sort_values([0, 1], ignore_index=True, inplace=True)
    tracks_df.loc[:, 3] = tracks_df.loc[:, 3].map(old_to_new_cm_map)

    content_df.loc[:, 0] = content_df.loc[:, 0].map(old_to_new_cm_map)
    content_df.loc[:, "index"] = content_df.loc[:, 0].copy()
    content_df.set_index("index", inplace=True)
    content_df.loc[tracks_df.loc[:, 3], :]

    tracks_df.to_csv(
        output_path + ".tracks.bed", sep="\t", index=False, header=False
    )
    content_df.loc[tracks_df.loc[:, 3], :].to_csv(
        output_path + ".content.txt", sep="\t", index=False, header=False
    )


