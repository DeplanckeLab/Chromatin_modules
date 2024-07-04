import itertools
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pyBigWig
import scipy


def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)


def get_corr_matrix(edgelist_df, output_path, dataset, gene, corr_file_name):
    df_sorted_pids = pd.DataFrame(
        index=list(
            set(edgelist_df.loc[:, "pid_1"]).union(set(edgelist_df.loc[:, "pid_2"]))
        )
    )
    df_sorted_pids.loc[:, "start"] = (
        df_sorted_pids.index.str.split(":").str[-2].astype(int)
    )
    df_sorted_pids.loc[:, "end"] = (
        df_sorted_pids.index.str.split(":").str[-1].astype(int)
    )
    df_sorted_pids = df_sorted_pids.sort_values(
        ["start", "end"], ascending=[True, True]
    )

    pivot_df = edgelist_df.pivot_table(
        index="pid_1", columns="pid_2", values="correlation", fill_value=np.nan
    )
    missing_col = list(set(pivot_df.index).difference(set(pivot_df.columns)))
    missing_idx = list(set(pivot_df.columns).difference(set(pivot_df.index)))

    pivot_df.loc[:, missing_col] = pd.DataFrame(
        np.zeros((pivot_df.shape[0], len(missing_col))), index=pivot_df.index
    )
    df_with_missing_rows = pd.DataFrame(
        np.zeros((len(missing_idx), pivot_df.shape[1])),
        index=missing_idx,
        columns=pivot_df.columns,
    )
    full_pivot_df = pd.concat([pivot_df, df_with_missing_rows], axis=0)
    del pivot_df
    full_pivot_df = full_pivot_df.loc[:, df_sorted_pids.index]
    full_pivot_df = full_pivot_df.loc[df_sorted_pids.index, :]
    # pivot_df = pivot_df.reindex(df_sorted_pids.index)
    np.fill_diagonal(
        full_pivot_df.values,
        pd.Series(data=[0] * full_pivot_df.shape[0], index=full_pivot_df.index),
    )
    for col in full_pivot_df.columns:
        missing_marks = full_pivot_df[full_pivot_df.loc[:, col].isnull()].index.tolist()
        if missing_marks:
            for row in missing_marks:
                full_pivot_df.loc[row, col] = full_pivot_df.loc[col, row]
    if not check_symmetric(full_pivot_df):
        print("Correlation matrix is not symmetric. Trying to fix it...")
        for r, row in enumerate(full_pivot_df.columns):
            for c, col in enumerate(list(full_pivot_df.columns)[r + 1 :]):
                if full_pivot_df.loc[row, col] != full_pivot_df.loc[col, row]:
                    if full_pivot_df.loc[row, col] == 0:
                        full_pivot_df.loc[row, col] = full_pivot_df.loc[col, row]
                    elif full_pivot_df.loc[col, row] == 0:
                        full_pivot_df.loc[col, row] = full_pivot_df.loc[row, col]
        print("Fixed!")
    if not os.path.exists(os.path.join(output_path, "corr_matrix_format", dataset)):
        os.makedirs(os.path.join(output_path, "corr_matrix_format", dataset))
    full_pivot_df.to_csv(
        os.path.join(
            output_path,
            "corr_matrix_format",
            dataset,
            "_".join([dataset, gene, corr_file_name]),
        ),
        sep="\t",
        index=True,
        header=True,
    )


def prepare_data_for_corr_heatmap(
    input_path, dataset, gene, corr_file_name, square_corr=True
):
    df_corr = pd.read_csv(
        os.path.join(
            input_path,
            "corr_matrix_format",
            dataset,
            "_".join([dataset, gene, corr_file_name]),
        ),
        sep="\t",
        index_col=0,
    )
    if square_corr:
        return np.square(df_corr)
    else:
        return df_corr


def prepare_data_for_hic_heatmap(input_path, dataset, gene, corr_file_name):
    hic_peak_df = pd.read_csv(
        os.path.join(
            input_path,
            "corr_matrix_format",
            dataset,
            "_".join([dataset, gene, corr_file_name]),
        ),
        sep="\t",
        index_col=0,
    )
    return np.log10(hic_peak_df + 1)


def get_cm_coordinates(corr_or_hic_peak_df, cm_dict, cms_to_plot, data_type="corr"):
    cm_coordinates = []
    for cm in cms_to_plot:
        if data_type == "corr":
            cm_borders = [
                corr_or_hic_peak_df.columns.get_loc(crd_peak)
                for crd_peak in cm_dict.get(cm)
            ]
        else:
            cm_borders = [
                corr_or_hic_peak_df.columns.get_loc(crd_peak.split("chr")[1])
                for crd_peak in cm_dict.get(cm)
            ]
        cm_coordinates.append([min(cm_borders), max(cm_borders)])
    return cm_coordinates


def get_cm_triangle_coordinates(pair, add_f=0.5):
    # The first side of a triangle
    x_1 = [
        -(pair[0] - add_f) * np.sqrt(2),
        -(
            (pair[0] - add_f) * np.sqrt(2)
            + ((pair[1] + add_f) - (pair[0] - add_f)) / np.sqrt(2)
        ),
    ]
    y_1 = [0, ((pair[1] + add_f) - (pair[0] - add_f)) / np.sqrt(2)]
    # Second side of a triangle
    x_2 = [
        -(
            (pair[0] - add_f) * np.sqrt(2)
            + ((pair[1] + add_f) - (pair[0] - add_f)) / np.sqrt(2)
        ),
        -(pair[1] + add_f) * np.sqrt(2),
    ]
    y_2 = [((pair[1] + add_f) - (pair[0] - add_f)) / np.sqrt(2), 0]
    return x_1, y_1, x_2, y_2


def plot_peak_corr_heatmap(
    dataset,
    corr_data,
    vmin_norm=0,
    vmax_norm=1,
    line_width=1,
    ylim_fraction=1,
    save_fig=False,
    square_corr=False,
    crd_coordinates=None,
    vcm_coordinates=None,
    phm_coordinates=None,
    cms_to_highlight=None,
    suffix=None,
    format_="pdf",
    cbar_label=None,
    output_path=None,
):
    if square_corr:
        cmap = plt.cm.Blues
    else:
        cmap = plt.cm.bwr
    fig, ax = plt.subplots(figsize=(20, 20))
    tr = matplotlib.transforms.Affine2D().rotate_deg(135)
    norm = matplotlib.colors.Normalize(vmin=vmin_norm, vmax=vmax_norm)
    add_f = 0.5
    if crd_coordinates != None:
        for pair in crd_coordinates:
            x_1, y_1, x_2, y_2 = get_cm_triangle_coordinates(pair, add_f=add_f)
            ax.plot(x_1, y_1, c="#A2A2A4", lw=line_width, ls="--")  # '#B0B0B2',
            ax.plot(x_2, y_2, c="#A2A2A4", lw=line_width, ls="--")  # '#B0B0B2',

    if vcm_coordinates != None:
        for pair in vcm_coordinates:
            x_1, y_1, x_2, y_2 = get_cm_triangle_coordinates(pair, add_f=add_f)
            ax.plot(x_1, y_1, c="#E9AD85", lw=line_width, ls="--")  # '#E8C6B0',
            ax.plot(x_2, y_2, c="#E9AD85", lw=line_width, ls="--")  # '#E8C6B0',

    if phm_coordinates != None:
        for pair in phm_coordinates:
            x_1, y_1, x_2, y_2 = get_cm_triangle_coordinates(pair, add_f=add_f)
            ax.plot(x_1, y_1, c="#68A19F", lw=line_width, ls="--")  # '#80B0AE',
            ax.plot(x_2, y_2, c="#68A19F", lw=line_width, ls="--")  # '#80B0AE',
    if cms_to_highlight != None:
        for pair in cms_to_highlight:
            x_1, y_1, x_2, y_2 = get_cm_triangle_coordinates(pair, add_f=add_f)
            ax.plot(x_1, y_1, c="darkorange", lw=line_width, ls="--")
            ax.plot(x_2, y_2, c="darkorange", lw=line_width, ls="--")

    ax.set_xlim((-corr_data.shape[0] * np.sqrt(2), 0))
    ax.set_ylim((0, corr_data.shape[0]))

    ax_pos = ax.get_position()
    fig_height = fig.get_size_inches()[1]
    ax_height = ax_pos.height * fig_height
    ax.set_ylim([0, ylim_fraction * ax_height])

    ax.set_xticks([])
    ax.set_yticks([])
    im = ax.imshow(
        corr_data,
        transform=tr + ax.transData,
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
    )
    ax.invert_xaxis()
    ax.axis("off")

    fig_colorbar, ax_colorbar = plt.subplots(1, 1, figsize=(1, 2))
    cbar = plt.colorbar(
        im,
        cax=ax_colorbar,
        orientation="vertical",
        norm=norm,  # define the colorbar range
    )
    ticks = np.linspace(vmin_norm, vmax_norm, 3)
    if vmax_norm <= 0.1:
        tick_labels = ["{:.3f}".format(x) for x in ticks]
    else:
        tick_labels = ["{:.1f}".format(x) for x in ticks]
    # Adjust cbar width, adapted from
    # https://stackoverflow.com/questions/53429258/matplotlib-change-colorbar-height-within-own-axes
    shrink = 0.1
    bbox_pos = ax_colorbar.get_position()
    new_width = bbox_pos.width * shrink
    pad = (bbox_pos.width - new_width) / 2
    bbox_pos.x0 = bbox_pos.x0 + pad
    bbox_pos.x1 = bbox_pos.x1 - pad
    ax_colorbar.set_position(bbox_pos)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(tick_labels, size=8)
    # cbar.ax.set_aspect(10)  # Set aspect ratio
    ax_colorbar.set_ylabel(cbar_label, size=10, rotation=270, labelpad=15)
    if save_fig:
        fig.savefig(
            os.path.join(
                output_path, "_".join([dataset, suffix, "peak_corr_heatmap." + format_])
            ),
            format=format_,
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )
        fig_colorbar.savefig(
            os.path.join(
                output_path,
                "_".join([dataset, suffix, "peak_corr_heatmap_colorbar." + format_]),
            ),
            format=format_,
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )


def plot_hic_peak_heatmap(
    resolution,
    dataset,
    hic_data,
    crd_coordinates,
    vcm_coordinates,
    phm_coordinates,
    line_width=1,
    vmax_norm=1,
    ylim_fraction=1,
    save_fig=False,
    suffix=None,
    cbar_label=None,
    output_path=None,
):
    cmap = plt.cm.Oranges
    fig, ax = plt.subplots(figsize=(20, 20))
    tr = matplotlib.transforms.Affine2D().rotate_deg(135)
    norm = matplotlib.colors.Normalize(vmin=0, vmax=vmax_norm)
    add_f = 0.5
    if crd_coordinates != None:
        for pair in crd_coordinates:
            x_1, y_1, x_2, y_2 = get_cm_triangle_coordinates(pair, add_f=add_f)
            ax.plot(x_1, y_1, c="#A2A2A4", lw=line_width, ls="--")  # '#B0B0B2',
            ax.plot(x_2, y_2, c="#A2A2A4", lw=line_width, ls="--")  # '#B0B0B2',

    if vcm_coordinates != None:
        for pair in vcm_coordinates:
            x_1, y_1, x_2, y_2 = get_cm_triangle_coordinates(pair, add_f=add_f)
            ax.plot(x_1, y_1, c="#E9AD85", lw=line_width, ls="--")  # '#E8C6B0',
            ax.plot(x_2, y_2, c="#E9AD85", lw=line_width, ls="--")  # '#E8C6B0',

    if phm_coordinates != None:
        for pair in phm_coordinates:
            x_1, y_1, x_2, y_2 = get_cm_triangle_coordinates(pair, add_f=add_f)
            ax.plot(x_1, y_1, c="#68A19F", lw=line_width, ls="--")  # '#80B0AE',
            ax.plot(x_2, y_2, c="#68A19F", lw=line_width, ls="--")  # '#80B0AE',

    ax.set_xlim((-hic_data.shape[0] * np.sqrt(2), 0))
    ax.set_ylim((0, hic_data.shape[0]))

    ax_pos = ax.get_position()
    fig_height = fig.get_size_inches()[1]
    ax_height = ax_pos.height * fig_height
    ax.set_ylim([0, ylim_fraction * ax_height])

    ax.set_xticks([])
    ax.set_yticks([])

    im = ax.imshow(
        hic_data,
        transform=tr + ax.transData,
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
    )

    ax.invert_xaxis()
    ax.axis("off")

    fig_colorbar, ax_colorbar = plt.subplots(1, 1, figsize=(1, 2))
    cbar = plt.colorbar(
        im,
        cax=ax_colorbar,
        orientation="vertical",
        norm=norm,  # define the colorbar range
    )
    ticks = np.linspace(0, vmax_norm, 3)
    if vmax_norm <= 0.1:
        tick_labels = ["{:.3f}".format(x) for x in ticks]
    else:
        tick_labels = ["{:.1f}".format(x) for x in ticks]
    # Adjust cbar width, adapted from
    # https://stackoverflow.com/questions/53429258/matplotlib-change-colorbar-height-within-own-axes
    shrink = 0.1
    bbox_pos = ax_colorbar.get_position()
    new_width = bbox_pos.width * shrink
    pad = (bbox_pos.width - new_width) / 2
    bbox_pos.x0 = bbox_pos.x0 + pad
    bbox_pos.x1 = bbox_pos.x1 - pad
    ax_colorbar.set_position(bbox_pos)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(tick_labels, size=8)
    # cbar.ax.set_aspect(10)  # Set aspect ratio
    ax_colorbar.set_ylabel(cbar_label, size=10, rotation=270, labelpad=15)

    if save_fig:
        fig.savefig(
            os.path.join(
                output_path,
                "_".join(
                    [
                        dataset,
                        suffix,
                        "peak_based_contact_matrix",
                        str(resolution),
                        "bp_interactions.pdf",
                    ]
                ),
            ),
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )
        fig_colorbar.savefig(
            os.path.join(
                output_path,
                "_".join(
                    [
                        dataset,
                        suffix,
                        "peak_based_contact_matrix",
                        str(resolution),
                        "bp_interactions_colorbar.pdf",
                    ]
                ),
            ),
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )


def plot_ruller_track(
    dataset,
    region_start,
    region_end,
    fig_width=20,
    save_fig=False,
    pid=None,
    output_path=None,
):
    fig, ax = plt.subplots(figsize=(fig_width, 1))
    # plot ruller track
    ax.hlines(0, region_start, region_end, color="black", lw=1)
    if fig_width < 15:
        if fig_width < 5:
            step = 2
        else:
            step = 3
    else:
        step = 6
    tick_list = np.linspace(region_start, region_end, step)
    for i in tick_list:
        if i == tick_list[0]:
            text_coord = round(i)
        elif i == tick_list[-1]:
            text_coord = round(
                region_end - 1 / (fig_width * 1.225) * (region_end - region_start)
            )
        else:
            text_coord = round(
                round(i) - 1 / (fig_width * 3) * (region_end - region_start)
            )
        ax.vlines(i, -0.15, 0.15, color="black", lw=1)
        # if i == region_end:
        str_ = "Mb"
        # else:
        #     str_ = ""
        ax.text(text_coord, 0.5, str(round(i / 10**6, 2)) + str_, fontsize=10)
    ax.axis("off")
    ax.set_ylim((-4, 0.5))
    ax.set_xlim((region_start, region_end))
    if save_fig:
        plt.savefig(
            os.path.join(output_path, "_".join([dataset, pid, "ruller_track.pdf"])),
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )


def plot_peaks(
    dataset,
    peaks_in_region,
    region_start,
    region_end,
    format_="string",
    peak_height=0.5,
    fig_width=20,
    save_fig=False,
    pid=None,
    output_path=None,
):
    fig, ax = plt.subplots(figsize=(fig_width, 1))
    if format_ == "string":
        peaks_in_region_ = [
            [
                peak_str.split(":")[0],
                int(peak_str.split(":")[-2]),
                int(peak_str.split(":")[-1]),
            ]
            for peak_str in peaks_in_region
        ]
    else:
        peaks_in_region_ = [["", peak_s, peak_e] for peak_s, peak_e in peaks_in_region]

    for i, (mark, peak_start, peak_end) in enumerate(peaks_in_region_):
        if format_ == "string":
            if mark == "H3K27ac":
                color = "#6D7A3E80"
            elif mark == "H3K4me1":
                color = "#4A4A9680"
        else:
            color = "black"
        rect = plt.Rectangle(
            (peak_start, 0 - peak_height / 2),
            peak_end - peak_start,
            peak_height,
            fill=True,
            color=color,
        )
        ax.add_patch(rect)

    ax.axis("off")
    ax.set_ylim((-4, 0.5))
    ax.set_xlim((region_start, region_end))
    if save_fig:
        plt.savefig(
            os.path.join(output_path, "_".join([dataset, pid, "peaks.pdf"])),
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )


def plot_genes(
    dataset,
    region_start,
    region_end,
    gene_exon_coordinates,
    ens_id_gene_symbol_mapping,
    gene_height=0.5,
    text_shift=-0.75,
    fig_width=20,
    save_fig=False,
    zoom_in=False,
    pid=None,
    output_path=None,
):
    fig, ax = plt.subplots(figsize=(fig_width, 1))
    shift = 0
    p_shift = 0
    for i, (ens_gene_id, (strand, exons)) in enumerate(gene_exon_coordinates.items()):
        exons = sorted(exons)
        gene_start = exons[0][0]
        gene_end = exons[-1][1]
        if i == 0:
            prev_gene_end = gene_end
        elif prev_gene_end < gene_start:
            p_shift += 1
        else:
            p_shift += 2

        gene_exon_starts, gene_exon_ends = [], []
        for exon_start, exon_end in exons:
            gene_exon_starts.append(exon_start)
            gene_exon_ends.append(exon_end)
            rect = plt.Rectangle(
                (exon_start, 0 - p_shift - gene_height / 2),
                exon_end - exon_start,
                gene_height,
                fill=True,
                color="black",
            )
            ax.add_patch(rect)
        ax.hlines(
            0 - p_shift,
            min(gene_exon_starts),
            max(gene_exon_ends),
            color="black",
            lw=1.5,
        )
        # if min(gene_exon_starts) < region_start:
        #     d = 1
        #     ax.plot((-d, +d), (-d - p_shift, d - p_shift), lw=1)  # top-left diagonal
        if strand == "+":
            marker = "->"
        else:
            marker = "-<"
        if zoom_in:
            if fig_width < 10:
                mark_every = 3
            else:
                mark_every = 2
        else:
            mark_every = 5
        gene_range = np.linspace(gene_start, gene_end, 20)
        ax.plot(
            gene_range[1:-1],
            [0 - p_shift] * len(gene_range[1:-1]),
            marker,
            lw=1.5,
            markevery=mark_every,
            markersize=5,
            color="black",
        )
        if len(ens_id_gene_symbol_mapping.get(ens_gene_id)) > 4:
            coeff = 0.06
        else:
            coeff = 0.045
        half_gene_name_len = (
            coeff
            * (region_end - region_start)
            / (fig_width * 1.225)
            * len(ens_id_gene_symbol_mapping.get(ens_gene_id))
        )
        if gene_start < region_start:
            x_text_coord = region_start
        elif gene_end > region_end:
            x_text_coord = region_end - 3 * half_gene_name_len
        else:
            x_text_coord = (gene_start + gene_end) / 2 - half_gene_name_len
        ax.text(
            x_text_coord,
            0 - p_shift + text_shift,
            ens_id_gene_symbol_mapping.get(ens_gene_id),
            fontsize=9.5,
            style="italic",
        )
        prev_gene_end = gene_end
        shift += 1 + p_shift
        if p_shift > 2:
            p_shift = 0
            shift = 0

    ax.axis("off")
    ax.set_ylim((-4, 0.5))
    ax.set_xlim((region_start, region_end))
    if save_fig:
        plt.savefig(
            os.path.join(output_path, "_".join([dataset, pid, "gene_tracks.pdf"])),
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )


def plot_CMs_rs(
    region_start,
    region_end,
    dataset,
    pid,
    save_fig=False,
    cm_height=0.35,
    fig_width=20,
    plot_variant=False,
    crd_color="#A2A2A4",
    vcm_color="#E9AD85",
    phm_color="#68A19F",
    crd_dict=None,
    vcm_dict=None,
    phm_dict=None,
    variant_location=None,
    variant_color="red",
    rs_ID=None,
    output_path=None,
    suffix=None,
):
    fig, ax = plt.subplots(figsize=(fig_width, 2))
    shift = 0
    if plot_variant:
        # plot rs variant
        if type(variant_location) in [int, float]:
            ax.vlines(variant_location, -0.5 - shift, -0.25 - shift, color="red", lw=2)
            ax.text(
                variant_location - len(rs_ID) * 100,
                -0.5 - 2 * shift,
                rs_ID,
                fontsize=10,
                color=variant_color,
            )
        else:
            for var_loc in variant_location:
                ax.vlines(var_loc, -0.75 - shift, -0.25 - shift, color="red", lw=2)
        shift += 5 * 0.5

    # plot CRDs
    if crd_dict is not None:
        for i, (_, crd_peaks) in enumerate(crd_dict.items()):
            crd_peaks = sorted(crd_peaks)
            crd_start = crd_peaks[0][0]
            crd_end = crd_peaks[-1][1]
            if i == 0:
                prev_crd_end = crd_end
            if prev_crd_end < crd_start:
                shift = prev_shift
            crd_peak_starts, crd_peak_ends = [], []
            for crd_peak_start, crd_peak_end in crd_peaks:
                crd_peak_starts.append(crd_peak_start)
                crd_peak_ends.append(crd_peak_end)
                rect = plt.Rectangle(
                    (crd_peak_start, 0 - shift - cm_height / 2),
                    crd_peak_end - crd_peak_start,
                    cm_height,
                    fill=True,
                    color=crd_color,
                )
                ax.add_patch(rect)

            ax.hlines(
                0 - shift,
                min(crd_peak_starts),
                max(crd_peak_ends),
                color=crd_color,
                lw=1.5,
            )
            prev_crd_end = crd_end
            prev_shift = shift
            shift += 0.5
    shift += 0.3
    # plot VCMs
    if vcm_dict is not None:
        for i, (_, vcm_peaks) in enumerate(vcm_dict.items()):
            vcm_peaks = sorted(vcm_peaks)
            vcm_start = vcm_peaks[0][0]
            vcm_end = vcm_peaks[-1][1]
            if i == 0:
                prev_vcm_end = vcm_end
            if prev_vcm_end < vcm_start:
                shift = prev_shift
            vcm_peak_starts, vcm_peak_ends = [], []
            for vcm_peak_start, vcm_peak_end in vcm_peaks:
                vcm_peak_starts.append(vcm_peak_start)
                vcm_peak_ends.append(vcm_peak_end)
                rect = plt.Rectangle(
                    (vcm_peak_start, 0 - shift - cm_height / 2),
                    vcm_peak_end - vcm_peak_start,
                    cm_height,
                    fill=True,
                    color=vcm_color,
                )
                ax.add_patch(rect)
            ax.hlines(
                0 - shift,
                min(vcm_peak_starts),
                max(vcm_peak_ends),
                color=vcm_color,
                lw=1.5,
            )
            prev_vcm_end = vcm_end
            prev_shift = shift
            shift += 0.5
    shift += 0.3
    # plot DAGs
    if phm_dict is not None:
        for i, (_, phm_peaks) in enumerate(phm_dict.items()):
            phm_peaks = sorted(phm_peaks)
            phm_start = phm_peaks[0][0]
            phm_end = phm_peaks[-1][1]
            if i == 0:
                prev_phm_end = phm_end
            if prev_phm_end < phm_start:
                shift = prev_shift
            phm_peak_starts, phm_peak_ends = [], []
            for phm_peak_start, phm_peak_end in phm_peaks:
                phm_peak_starts.append(phm_peak_start)
                phm_peak_ends.append(phm_peak_end)
                rect = plt.Rectangle(
                    (phm_peak_start, 0 - shift - cm_height / 2),
                    phm_peak_end - phm_peak_start,
                    cm_height,
                    fill=True,
                    color=phm_color,
                )
                ax.add_patch(rect)
            ax.hlines(
                0 - shift,
                min(phm_peak_starts),
                max(phm_peak_ends),
                color=phm_color,
                lw=1.5,
            )
            prev_phm_end = phm_end
            prev_shift = shift
            shift += 0.5

    ax.axis("off")
    ax.set_ylim((-7, 0.5))
    ax.set_xlim((region_start, region_end))
    if suffix is None:
        suffix = "clomics_vcmtools_phm"
    if save_fig:
        plt.savefig(
            os.path.join(output_path, "_".join([dataset, pid, suffix, "tracks.pdf"])),
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )


def plot_TADs(
    dataset,
    region_start,
    region_end,
    tad_list,
    fig_width=20,
    save_fig=False,
    pid=None,
    output_path=None,
):
    fig, ax = plt.subplots(figsize=(fig_width, 1))
    shift = 0
    prev_tad_end = tad_list[0][1]
    for i, (tad_start, tad_end) in enumerate(tad_list):
        if prev_tad_end < tad_start:
            shift = prev_shift
        ax.hlines(3 - shift, tad_start, tad_end, color="#B0B0B0", lw=6)
        prev_tad_end = tad_end
        prev_shift = shift
        shift += 0.4
    ax.axis("off")
    ax.set_xlim((region_start, region_end))
    ax.set_ylim((-1, 5))
    if save_fig:
        plt.savefig(
            os.path.join(output_path, "_".join([dataset, pid, "TAD_tracks.pdf"])),
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )


def plot_CTCF(
    dataset,
    region_start,
    region_end,
    forward_ctcf_positions,
    reverse_ctcf_positions,
    fig_width=20,
    save_fig=False,
    output_path=None,
    pid=None,
):
    fig, ax = plt.subplots(figsize=(fig_width, 2))
    # plot forward CTCF
    fw_x = (
        [region_start - 500]
        + [fw[0] for fw in forward_ctcf_positions]
        + [region_end + 500]
    )
    ax.plot(
        fw_x[1:-1],
        [0.75] * len(fw_x[1:-1]),
        "->",
        lw=0.8,
        markersize=8,
        color="#225F71",
    )
    ax.plot(
        fw_x,
        [0.75] * len(fw_x),
        lw=0.8,
        color="#225F71",
    )
    # plot reverse CTCF
    rv_x = (
        [region_start - 500]
        + [rv[0] for rv in reverse_ctcf_positions]
        + [region_end + 500]
    )
    ax.plot(
        rv_x[1:-1],
        [0.5] * len(rv_x[1:-1]),
        "-<",
        lw=0.8,
        markersize=8,
        color="#EA3449",
    )
    ax.plot(
        rv_x,
        [0.5] * len(rv_x),
        lw=0.8,
        color="#EA3449",
    )
    ax.axis("off")
    ax.set_xlim((region_start, region_end))
    ax.set_ylim((0, 1.5))
    if save_fig:
        plt.savefig(
            os.path.join(
                output_path, "_".join([dataset, pid, "CTCF_forward_reverse_tracks.pdf"])
            ),
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )


def get_binirazed_bw(bw_file, chromosome, start, end, n_bins):
    bw_values = bw_file.values(chromosome, start, end)
    return [
        np.mean(bw_values[i : i + n_bins]) for i in range(0, len(bw_values), n_bins)
    ]


def get_max_signal_bigwig(chromosome, region_start, region_end, bw_path, file_name):
    max_signal = 0
    # try:
    bw = pyBigWig.open(os.path.join(bw_path, file_name))
    # except:
    # return None
    n = int(np.ceil(0.05 / 24.5 * (region_end - region_start)))
    bw_bins = get_binirazed_bw(bw, chromosome, region_start, region_end, n_bins=n)
    bw.close()

    max_signal = max(max(bw_bins), max_signal)
    return max_signal


def plot_several_bigwig_profiles(
    dataset,
    chromosome,
    region_start,
    region_end,
    sample_dict,
    n_samples,
    bw_path,
    k27ac_file_name,
    k4me1_file_name,
    fig_width=20,
    save_fig=False,
    output_path=None,
    pid=None,
):
    if fig_width < 10:
        coeff = 1.3
    else:
        coeff = 2
    fig, ax = plt.subplots(n_samples * 2, 1, figsize=(fig_width, n_samples * coeff))
    if fig_width < 15:
        step = 3
    else:
        step = 6
    max_k27ac_signal = 0
    max_k4me1_signal = 0
    i = 0
    for _, sample_list in sample_dict.items():
        counter = 0
        for sample_id in sample_list:
            try:
                k27ac_bw = pyBigWig.open(
                    os.path.join(bw_path, "_".join([sample_id, k27ac_file_name]))
                )
                k4me1_bw = pyBigWig.open(
                    os.path.join(bw_path, "_".join([sample_id, k4me1_file_name]))
                )
                counter += 1
            except:
                continue
            if counter == 3:
                continue
            if fig_width < 10:
                n = int(np.ceil(0.08 / 24.5 * (region_end - region_start)))
            else:
                n = int(np.ceil(0.05 / 24.5 * (region_end - region_start)))

            k27ac_bw_bins = get_binirazed_bw(
                k27ac_bw, chromosome, region_start, region_end, n_bins=n
            )
            k4me1_bw_bins = get_binirazed_bw(
                k4me1_bw, chromosome, region_start, region_end, n_bins=n
            )
            k27ac_bw.close()
            k4me1_bw.close()

            max_k27ac_signal = max(max(k27ac_bw_bins), max_k27ac_signal)
            max_k4me1_signal = max(max(k4me1_bw_bins), max_k4me1_signal)
            ax[i].bar(
                np.arange(len(k27ac_bw_bins)),
                k27ac_bw_bins,
                color="#6D7A3E",
                rasterized=True,
                lw=0,
            )
            ax[i + 1].bar(
                np.arange(len(k4me1_bw_bins)),
                k4me1_bw_bins,
                color="#4A4A96",
                rasterized=True,
                lw=0,
            )
            ax[i].set_ylabel(sample_id)
            # ax[i + 1].set_ylabel(sample_id)
            i += 2
    for i in range(0, n_samples * 2):
        if i % 2 == 0:
            # ax[i].set_ylim((0, max_k27ac_signal + 1))
            ax[i].set_ylim((0, max_k27ac_signal - 0.1 * max_k27ac_signal))
            ax[i].set_xlim((0, len(k27ac_bw_bins)))
        else:
            # ax[i].set_ylim((0, max_k4me1_signal + 1))
            ax[i].set_ylim((0, max_k4me1_signal - 0.1 * max_k4me1_signal))
            ax[i].set_xlim((0, len(k4me1_bw_bins)))
        ax[i].spines["right"].set_visible(False)
        ax[i].spines["top"].set_visible(False)

        ax[i].set_xticks(
            np.linspace(0, len(k27ac_bw_bins), step),
            ["" for _ in np.linspace(region_start, region_end, step)],
        )
    label_patches = [
        matplotlib.patches.Patch(color="#6D7A3E", alpha=1, label="H3K27ac"),
        matplotlib.patches.Patch(color="#4A4A96", alpha=1, label="H3K4me1"),
    ]
    # box = ax[-1].get_position()
    # ax[-1].set_position([
    #     box.x0,
    #     box.y0,
    #     box.width,
    #     box.height
    # ])
    ax[-1].legend(
        handles=label_patches,
        loc="lower center",
        bbox_to_anchor=(0.5, -1),
        fancybox=False,
        shadow=False,
        ncol=2,
    )
    fig.tight_layout(h_pad=0.5)
    if save_fig:
        plt.savefig(
            os.path.join(
                output_path, "_".join([dataset, pid, "K27ac_k4me1_bigwig_tracks.pdf"])
            ),
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )


def plot_bigwig_profile(
    dataset,
    chromosome,
    region_start,
    region_end,
    sample_id,
    file_name,
    mark,
    bw_path,
    fig_width=20,
    bin_factor=0.05,
    save_fig=False,
    output_path=None,
    pid=None,
    y_max=None,
    is_left=True,
    color=None,
):
    max_signal = 0
    max_signal = 0
    try:
        bw = pyBigWig.open(os.path.join(bw_path, file_name))
    except:
        return None
    n = int(np.ceil(bin_factor / 24.5 * (region_end - region_start)))
    bw_bins = get_binirazed_bw(bw, chromosome, region_start, region_end, n_bins=n)
    bw.close()
    max_signal = max(max(bw_bins), max_signal)

    fig, ax = plt.subplots(1, 1, figsize=(fig_width, 0.5))
    if fig_width < 15:
        step = 3
    else:
        step = 6
    if color is None:
        if mark == "H3K27ac":
            color = "#6D7A3E"
        elif mark == "H3K4me1":
            color = "#4A4A96"
    ax.bar(np.arange(len(bw_bins)), bw_bins, color=color, rasterized=True, lw=0)
    # if is_left:
    ax.set_ylabel(sample_id)
    if y_max is None:
        ax.set_ylim((0, max_signal + 0.1 * max_signal))
    else:
        if y_max < 5:
            y_max = 5
        ax.set_ylim((0, y_max))
    ax.set_xlim((0, len(bw_bins)))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    if not is_left:
        ax.spines["left"].set_visible(False)
        ax.set_yticks([])
        ax.set_yticklabels("")
    ax.set_xticks(
        np.linspace(0, len(bw_bins), step),
        ["" for _ in np.linspace(region_start, region_end, step)],
    )
    #     ax.legend(
    #         handles=label_patches,
    #         loc="lower center",
    #         bbox_to_anchor=(0.5, -1),
    #         fancybox=False,
    #         shadow=False,
    #         ncol=1,
    #     )
    if save_fig:
        plt.savefig(
            os.path.join(
                output_path,
                "_".join([dataset, pid, sample_id, mark, "bigwig_tracks.pdf"]),
            ),
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )


def get_bp_profile_from_bigwig(
    chromosome,
    region_start,
    region_end,
    file_name_list,
    bw_path,
    ax,
    color=None,
    method="CV",
):
    n = int(np.ceil(0.05 / 24.5 * (region_end - region_start)))

    bw_profiles_per_sample = []
    for file_name in file_name_list:
        try:
            bw = pyBigWig.open(os.path.join(bw_path, file_name))
            bw_bins = get_binirazed_bw(
                bw, chromosome, region_start, region_end, n_bins=n
            )
            bw_profiles_per_sample.append(np.array(bw_bins))
            bw.close()
        except:
            continue

    std_per_bin = np.std(bw_profiles_per_sample, axis=0)
    if method == "CV":
        average_per_bin = np.mean(bw_profiles_per_sample, axis=0)
        cv_per_bin = std_per_bin / average_per_bin
        cv_per_bin[cv_per_bin == np.inf] = 0
        y_label = "CV"
    else:
        cv_per_bin = std_per_bin
        y_label = "$\\sigma$"

    if color is None:
        color = "gray"
    ax.plot(
        np.arange(len(cv_per_bin)),
        cv_per_bin,
        color=color,
        rasterized=True,
    )
    ax.set_ylabel("Bin-based " + y_label, size=12)
    ax.set_xlim((0, len(bw_bins)))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xticks(
        np.linspace(0, len(bw_bins), 6),
        ["" for _ in np.linspace(region_start, region_end, 6)],
    )


def plot_arcs(
    dataset,
    chromosome,
    corr_or_hic_peak_df,
    common_peaks,
    region_start,
    region_end,
    interaction_threshold=0.2,
    save_fig=False,
    fig_width=20,
    cmap=plt.cm.Reds,
    min_val=-1,
    max_val=1,
    output_path=None,
    pid=None,
):
    cmap = cmap
    min_corr = min_val
    max_corr = max_val
    norm = matplotlib.colors.TwoSlopeNorm(vmin=min_corr, vcenter=0, vmax=max_corr)
    if fig_width < 10:
        fig_height = 1.3
    else:
        fig_height = 2
    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
    corr_or_hic_peak_df.index = corr_or_hic_peak_df.index.str.split(
        chromosome + ":"
    ).str[1]
    corr_or_hic_peak_df.columns = corr_or_hic_peak_df.columns.str.split(
        chromosome + ":"
    ).str[1]
    for peak_1, peak_2 in itertools.combinations(sorted(common_peaks), 2):
        try:
            peak_1_start, peak_1_end = peak_1
            peak_2_start, peak_2_end = peak_2
            peak_1_center = (peak_1_start + peak_1_end) / 2
            peak_2_center = (peak_2_start + peak_2_end) / 2
            peak_correlation = corr_or_hic_peak_df.loc[
                str(peak_1_start) + ":" + str(peak_1_end),
                str(peak_2_start) + ":" + str(peak_2_end),
            ]
            if abs(peak_correlation) < interaction_threshold:
                continue
            if peak_1_center <= peak_2_center:
                left_center = peak_1_center
                right_center = peak_2_center
            else:
                left_center = peak_2_center
                right_center = peak_1_center
            x = [left_center, (left_center + right_center) / 2, right_center]
            y = [0, (right_center - left_center) / (10**2), 0]
            support = np.linspace(left_center, right_center, 1000)
            spl = scipy.interpolate.CubicSpline(x, y)
            y_smooth = spl(support)
            plt.plot(
                support,
                -y_smooth,
                linewidth=1,  # peak_correlation * np.sqrt(3.5),
                color=cmap(norm(peak_correlation)),
                alpha=0.5,
            )
        except:
            continue
    # ax.hlines(0, region_start, region_end, color="black", lw=1)
    ax.set_xlim((region_start, region_end))
    # plt.ylim((-9, 0))
    ax.axis("off")
    if save_fig:
        plt.savefig(
            os.path.join(output_path, "_".join([dataset, pid, "CM_arc_plots.pdf"])),
            bbox_inches="tight",
            edgecolor=None,
            facecolor=None,
            transparent=True,
            dpi=300,
        )
