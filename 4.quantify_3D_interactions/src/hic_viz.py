import collections
import cooler
import matplotlib
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import numpy as np
import pandas as pd
import pylab as pl
import scipy.stats as stats
import seaborn as sns
import statsmodels.stats.multitest as multi
import sys

# import straw
import tqdm
import warnings

from os import listdir
from os.path import join
from matplotlib import cm, transforms
from scipy.sparse import coo_matrix
from scipy.stats import mannwhitneyu
from scipy import stats
from scipy.interpolate import CubicSpline
from sklearn.decomposition import PCA


class HiVIZ:
    def __init__(self, numpy_mtx):
        self.numpy_mtx = self.numpy_mtx

    def get_rotation_coord(self):
        n = self.numpy_mtx.shape[1]
        d = self.numpy_mtx.copy()
        mask = np.zeros_like(d, dtype=np.bool)
        mask[np.tril_indices_from(mask)] = True
        C = np.ma.masked_array(d, mask=mask)

        # Transformation matrix for rotating the heatmap.
        xy = np.array([(x, y) for x in range(n + 1) for y in range(n + 1)])
        rot = np.array([[1, -1], [1, 1]]) * np.sqrt(2) / 2
        A = np.dot(xy, rot)

        # Plot the correlation heatmap triangle.
        X = A[:, 1].reshape(n + 1, n + 1)
        Y = A[:, 0].reshape(n + 1, n + 1)
        return X, Y, C

    def gene_coordinates(
        self,
        region_start,
        region_end,
        new_coord=False,
        gene_dict=None,
        old_gene_dict=None,
        old_region_end=None,
    ):
        plt.figure(figsize=(15, 1))
        plt.hlines(0, region_start, region_end, color="black", lw=1)
        for i in np.arange(region_start, region_end + 1000, 1000):
            if (i != region_start) and (i != region_end) and (i % (1 * 10**4) == 0):
                plt.vlines(
                    i, -(10 ** (-7.3)) / 4, 10 ** (-7.3) / 4, color="black", lw=1
                )
                plt.text(
                    i - 1200,
                    -(10 ** (-7.3)),
                    str(round(i / 10**6, 2)) + " Mb",
                    fontsize=10,
                )
            if i == region_start:
                plt.vlines(
                    i, -(10 ** (-7.3)) / 4, 10 ** (-7.3) / 4, color="black", lw=1
                )
                plt.text(
                    i - 1200,
                    -(10 ** (-7.3)),
                    str(round(i / 10**6, 2)) + " Mb",
                    fontsize=10,
                )

            elif i == region_end:
                plt.vlines(
                    i, -(10 ** (-7.3)) / 4, 10 ** (-7.3) / 4, color="black", lw=1
                )
                if new_coord:
                    plt.text(
                        i - 1200,
                        -(10 ** (-7.3)),
                        str(round(old_region_end / 10**6, 2)) + " Mb",
                        fontsize=10,
                    )
                else:
                    plt.text(
                        i - 1200,
                        -(10 ** (-7.3)),
                        str(round(i / 10**6, 2)) + " Mb",
                        fontsize=10,
                    )
        for gene, gene_coordinates in gene_dict.items():
            if gene == "DAG":
                plt.hlines(
                    0,
                    gene_coordinates[0] + np.sqrt(2) / 2,
                    gene_coordinates[1] + np.sqrt(2) / 2,
                    color="lightblue",
                    lw=5,
                )
                plt.text(
                    gene_coordinates[0]
                    + (gene_coordinates[1] - gene_coordinates[0]) / 2
                    + np.sqrt(2) / 2,
                    10 ** (-7.5),
                    gene,
                    fontsize=10,
                )
            elif gene == "Lead":
                plt.hlines(
                    0,
                    gene_coordinates[0] + np.sqrt(2) / 2,
                    gene_coordinates[1] + np.sqrt(2) / 2,
                    color="red",
                    lw=5,
                )
                plt.text(gene_coordinates[0], 10 ** (-7.5), gene, fontsize=10)
            elif gene == "Dependent":
                for dependent_peak in gene_coordinates:
                    plt.hlines(
                        0,
                        dependent_peak[0] + np.sqrt(2) / 2,
                        dependent_peak[1] + np.sqrt(2) / 2,
                        color="navy",
                        lw=5,
                    )
        plt.xlim((region_start, region_end))
        plt.ylim((-(10 ** (-7)), 10 ** (-7)))
        plt.axis("off")
        if zoom_in:
            zoom_type = "zoom_in"
        else:
            zoom_type = "zoom_out"
        if new_coord:
            suffix = "_without_zeros"
        else:
            suffix = "_with_zeros"

    #     output_path = '/home/pushkare/SVFASRAW/pushkare/filtered_mm10_Adcy3_enhancer_interactions/plots/log10_norm_reads'
    #     plt.savefig(join(output_path, zoom_type, 'gene_coordinates_mm10_Prkarb2_' + zoom_type + suffix + '.pdf'),
    #                 bbox_inches='tight',
    #                 edgecolor=None,
    #                 facecolor=None,
    #                 dpi=500)

    def plot_arcs(
        self, region_start, region_end, d_scale, quantile=False, new_coord=False, q=None
    ):
        if quantile:
            data_quantile = np.quantile(self.numpy_mtx, q)
        row_idx = np.arange(region_start, region_end, 1000) / 1000
        row_idx = (row_idx - np.min(row_idx)).astype(int)
        region = np.linspace(region_start, region_end, self.numpy_mtx.shape[0] - 1)
        dict_with_idx = dict(list(zip(region, row_idx)))

        plt.figure(figsize=(15, 10))
        plt.hlines(0, region_start, region_end, color="black", lw=1)
        for i, el_i in enumerate(region[:-1]):
            for el_j in region[i + 1 :]:
                d = self.numpy_mtx[dict_with_idx[el_i], dict_with_idx[el_j]]
                if quantile:
                    if d > data_quantile:
                        x = [el_i, (el_i + el_j) / 2, el_j]
                        y = [
                            0,
                            (el_j - el_i) / (20 * (10**2.5)),
                            0,
                        ]  # 65), 0] # 4.6), 0], zoom in
                        support = np.linspace(el_i, el_j, 10000)
                        spl = CubicSpline(x, y)
                        y_smooth = spl(support)
                        plt.plot(
                            support,
                            -y_smooth,
                            linewidth=d**2,  # for log / 60, #
                            color=cm.gist_heat_r(d**2 / d_scale),
                        )
                    else:
                        continue
                else:
                    x = [el_i, (el_i + el_j) / 2, el_j]
                    y = [0, (el_j - el_i) / (10**4.2), 0]  # 65), 0] # 4.3), 0], zoom in
                    support = np.linspace(el_i, el_j, 10000)
                    spl = CubicSpline(x, y)
                    y_smooth = spl(support)
                    plt.plot(
                        support,
                        -y_smooth,
                        linewidth=d**2,  # for log  / 60, #
                        color=cm.gist_heat_r(d**2 / d_scale),
                    )  # , for log  #/ 600)) #
        for j in region:
            plt.scatter(j, 0, marker="o", s=50, color="black")
        plt.xlim((region_start, region_end))
        plt.ylim((-10, 0))
        plt.axis("off")

    #     output_path = '/home/pushkare/SVFASRAW/pushkare/filtered_mm10_Adcy3_enhancer_interactions/plots/log10_norm_reads'
    #     plt.savefig(join(output_path, zoom_type, 'arcs_mm10_Prkarb2_' + zoom_type + suffix1 + suffix2 + '.pdf'),
    #                 bbox_inches='tight',
    #                 edgecolor=None,
    #                 facecolor=None,
    #                 dpi=500)
