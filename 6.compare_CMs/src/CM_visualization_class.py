import itertools
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import networkx as nx
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys
import warnings

warnings.filterwarnings("ignore")

class CM_VIZ:   

    def __init__(
            self,
            input_path=None,
            methods_list=None,
            score_type=None,
            similarity_threshold=None,
            pp_threshold=None,
            use_cmQTLs=False,
            save_figure=False
        ):
        self.input_path = input_path
        self.methods_list = methods_list
        self.score_type = score_type
        self.use_cmQTLs = use_cmQTLs
        if self.use_cmQTLs:
            self.folder_name = 'modules_with_QTL'
        else:
            self.folder_name = 'all_modules'
        self.save_figure = save_figure
        self.similarity_threshold = similarity_threshold
        self.pp_threshold = pp_threshold

    def load_npy(
            self,
            method_A,
            method_B,
            dataset,
            overlap_dict_path
        ):
        prefix = '_'.join([method_A, method_B])
        if self.use_cmQTLs:
            suffix = 'scores_dict_only_with_QTL.npy'
        else:
            suffix = 'scores_dict.npy'
        return np.load(
            os.path.join(
                overlap_dict_path,
                '_'.join([prefix, dataset, self.score_type, suffix])
            ),
            allow_pickle=True
        ).item()

    def load_tracks(
            self,
            method,
            dataset,
            similarity_threshold,
            tracks_output_path
        ):
        if (similarity_threshold > 0) and (similarity_threshold < 1):
            suffix = 'greater'
        elif (similarity_threshold == 0) or (similarity_threshold == 1):
            suffix = 'equal'
        elif similarity_threshold is None:
            suffix = 'all'
        similarity_threshold = np.round(similarity_threshold, 3)
        return pd.read_csv(
            os.path.join(
                tracks_output_path,
                '_'.join([
                    dataset,
                    method,
                    self.score_type,
                    'scores',
                    suffix,
                    str(similarity_threshold) + '.tracks.bed'])
                ),
            sep='\t',
            header=None
            )

    def get_CM_overlap_df(
            self,
            dataset,
            threshold,
            tracks_output_path,
            tracks_paths_dict,
            pairwise
        ):
        universal_CM_dict = {}
        for method in self.methods_list:
            universal_CM_dict[method] = self.load_tracks(
                    method,
                    dataset,
                    threshold,
                    tracks_output_path
                    )
        if pairwise:
            method_pairs = list(itertools.combinations(self.methods_list, 2)) 
            jaccard_list = []
            for method_A, method_B in method_pairs:
                n_CMs_A = pd.read_csv(tracks_paths_dict[method_A], sep='\t', header=None).shape[0]
                n_CMs_B = pd.read_csv(tracks_paths_dict[method_B], sep='\t', header=None).shape[0]
                n_CMs = n_CMs_A + n_CMs_B
                overlap_size = universal_CM_dict[method_A].shape[0] + universal_CM_dict[method_B].shape[0]
                jaccard_list.append([method_A, method_B, overlap_size / (n_CMs - overlap_size)])
                jaccard_list.append([method_B, method_A, overlap_size / (n_CMs - overlap_size)])
            jaccard_df = pd.DataFrame(jaccard_list, columns=['from', 'to', 'weight'])
            G = nx.Graph()
            G.add_nodes_from(set(jaccard_df['from']).union(set(jaccard_df['to'])))
            for node_1, node_2, weight in zip(jaccard_df['from'], jaccard_df['to'], jaccard_df['weight']):
                G.add_edge(node_1, node_2, weight=weight)
            return nx.to_pandas_adjacency(G, dtype=float)
        else:
            n_CMs = np.sum([
                pd.read_csv(tracks_path, sep='\t', header=None).shape[0]
                for _, tracks_path in tracks_paths_dict.items()
            ])
            overlap_size = np.sum([
                len(universal_CM_list)
                for _, universal_CM_list in universal_CM_dict.items()
            ])
            return pd.DataFrame.from_dict(
                {
                    '_'.join(list(tracks_paths_dict.keys())): overlap_size / (n_CMs - overlap_size)
                    },
                orient='index')
    
    def cm_based_jaccard_df_for_similarity_threshold(
            self,
            dataset,
            tracks_paths_dict,
            pairwise=False
        ):
        if self.pp_threshold is not None:
            str_ = '/pp_threshold_' + str(self.pp_threshold)
        else:
            str_ = ''
        tracks_output_path = os.path.join(
            self.input_path,
            dataset,
            self.folder_name + str_,
            self.score_type,
            'subset_tracks_and_content',
            self.folder_name
        )
        if type(self.similarity_threshold) == float:
            self.similarity_threshold = np.round(self.similarity_threshold, 3)
            return self.get_CM_overlap_df(
                dataset,
                self.similarity_threshold,
                tracks_output_path,
                tracks_paths_dict,
                pairwise
            )
        else:
            dfs_per_threshod = {}
            for threshold in self.similarity_threshold:
                threshold = np.round(threshold, 3)
                dfs_per_threshod[threshold] = self.get_CM_overlap_df(
                    dataset,
                    threshold,
                    tracks_output_path,
                    tracks_paths_dict,
                    pairwise
                )
            return dfs_per_threshod
        
    def eCDF_plot(
            self,
            dataset,
            cmap,
            title=None
        ):
        all_method_pairs = list(itertools.combinations(self.methods_list, 2))
        if self.pp_threshold is not None:
            str_ = 'pp_threshold_' + str(self.pp_threshold)
        else:
            str_ = ''
        plt.figure(figsize=(7, 7))
        for i, (method_A, method_B) in enumerate(all_method_pairs):
            A_B_F1_scores = self.load_npy(
                method_A,
                method_B,
                dataset,
                os.path.join(
                    self.input_path,
                    dataset,
                    self.score_type,
                    self.score_type + '_scores'
                )
            )
            sim_scores_df = pd.DataFrame(
                [
                    [method_A + '_' + cm_ids_A_B.split('_')[0],
                    method_B + '_' + cm_ids_A_B.split('_')[1],
                    sim_score]
                    for cm_ids_A_B, sim_score in A_B_F1_scores.items()
                ],
                columns=['cms_A', 'cms_B', 'score']
            )
            sim_scores_df = sim_scores_df[~sim_scores_df.loc[:, 'cms_A'].str.contains('NAN')]
            sim_scores_df = sim_scores_df[~sim_scores_df.loc[:, 'cms_B'].str.contains('NAN')]
            sns.kdeplot(
                list(sim_scores_df.loc[:, 'score']),
                cumulative=True,
                label=method_A + ' vs ' + method_B,
                color=cmap[i]
            )
            plt.axvline(
                np.mean(list(sim_scores_df.loc[:, 'score'])),
                linestyle='--',
                color=cmap[i],
                linewidth=1
            )

        handles, labels = plt.gca().get_legend_handles_labels()
        empty_patch = patches.Patch(
                linestyle='--',
                edgecolor='grey',
                linewidth=1,
                fill=False,
                label='Average score'
            )  # create a patch with no color

        handles.append(empty_patch)  # add new patches and labels to list
        labels.append('Average score')
        plt.legend(handles, labels)
        plt.xlim((0, 1))
        plt.ylim((0, 1))
        plt.xlabel('Reproducibility score', size=12)
        plt.ylabel('Density', size=12)
        if self.use_cmQTLs:
            if title is None:
                plt.title(
                    'Peak-based cell-type-specificity scores\nfor chromatin modules with cmQTL in ' + dataset,
                    size=12
                )
            else:
                plt.title(title, size=12)
        else:
            if title is None:
                plt.title(
                    'Peak-based cell-type-specificity\nfor chromatin modules in ' + dataset,
                    size=12
                )
            else:
                plt.title(title, size=12)

        if self.save_figure:
            if self.pp_threshold is not None:
                str_ = 'pp_threshold_' + str(self.pp_threshold)
            else:
                str_ = ''
            figure_out_path = os.path.join(
                self.input_path,
                'plots',
                dataset,
                self.folder_name,
                self.score_type + str_,
            )
            if not os.path.exists(figure_out_path):
                os.makedirs(figure_out_path)
            if self.use_cmQTLs:
                plt.savefig(
                    os.path.join(
                        figure_out_path,
                        '_'.join([dataset] + self.methods_list  + ['with_cmQTLs_eCDF.pdf'])
                        ),
                    dpi=300,
                    transparent=True,
                    bbox_inches='tight'
                )
            else:
                plt.savefig(
                    os.path.join(
                        figure_out_path,
                        '_'.join([dataset] + self.methods_list  + ['all_CMs_eCDF.pdf'])
                        ),
                    dpi=300,
                    transparent=True,
                    bbox_inches='tight'
                )

    def hex_plot(
        self,
        dataset
    ):
        all_method_pairs = list(itertools.permutations(self.methods_list, 2))
        if self.pp_threshold is not None:
            str_ = 'pp_threshold_' + str(self.pp_threshold)
        else:
            str_ = ''
        scores_list = []
        for i, (method_A, method_B) in enumerate(all_method_pairs):
            A_B_F1_scores = self.load_npy(
                method_A,
                method_B,
                dataset,
                os.path.join(
                    self.input_path,
                    dataset,
                    self.score_type,
                    self.score_type + '_scores'
                )
            )
            scores_list.append(list(A_B_F1_scores.values()))
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        y = 10
        ax.hexbin(
            scores_list[0],
            scores_list[1],
            gridsize=(int(np.sqrt(3) * y), y),
            cmap='Blues',
            bins='log'
        )
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.axvline(
            np.mean(
                scores_list[0]
            ),
            linestyle='--',
            color='black',
            linewidth=1
        )
        ax.axvline(
            np.mean(
                [el for el in scores_list[0] if el > 0]
            ),
            linestyle='--',
            color='gray',
            linewidth=1
        )
        ax.axhline(
            np.mean(
                scores_list[0]
            ),
            linestyle='--',
            color='black',
            linewidth=1
        )
        ax.axhline(
            np.mean(
                [el for el in scores_list[0] if el > 0]
            ),
            linestyle='--',
            color='gray',
            linewidth=1
        )

    def grid_plot(
        self,
        dataset,
        method_A,
        method_B,
        method_color_dict 
    ):
        A_B_F1_scores = self.load_npy(
            method_A,
            method_B,
            dataset,
            os.path.join(
                self.input_path,
                dataset,
                self.score_type,
                self.score_type + '_scores'
            )
        )
        sim_scores_df = pd.DataFrame(
            [
                [method_A + '_' + cm_ids_A_B.split('_')[0],
                method_B + '_' + cm_ids_A_B.split('_')[1],
                sim_score]
                for cm_ids_A_B, sim_score in A_B_F1_scores.items()
            ],
            columns=['cms_A', 'cms_B', 'score']
        )
        sim_scores_df = sim_scores_df[~sim_scores_df.loc[:, 'cms_A'].str.contains('NAN')]
        sim_scores_df = sim_scores_df[~sim_scores_df.loc[:, 'cms_B'].str.contains('NAN')]
        # https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html
        fig = plt.figure(figsize=(6, 6))
        gs = fig.add_gridspec(
            2,
            2,
            width_ratios=(3, 1),
            height_ratios=(1, 3),
            wspace=0.05,
            hspace=0.05
        )

        # Create the Axes.
        ax = fig.add_subplot(gs[1, 0])
        ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
        ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

        # Draw the scatter plot and marginals.
        # no labels
        ax_histx.tick_params(axis="x", labelbottom=False)
        ax_histy.tick_params(axis="y", labelleft=False)

        x = list(sim_scores_df.loc[:, 'score'])
        y = list(sim_scores_df.loc[:, 'score'])
        # the scatter plot:
        ax.scatter(
            x,
            y,
            c='lightgray',
            alpha=0.5,
            edgecolors='gray',
            linewidths=1
        )
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)

        # now determine nice limits by hand:
        binwidth = 0.05
        bins = np.arange(0, 1 + binwidth, binwidth)
        ax_histx.hist(
            x,
            bins=bins,
            density=True,
            color=method_color_dict[method_A]
        )
        ax_histy.hist(
            y,
            bins=bins,
            orientation='horizontal',
            density=True,
            color=method_color_dict[method_B]
        )
        ax_histx.spines.right.set_visible(False)
        ax_histx.spines.top.set_visible(False)
        ax_histy.spines.right.set_visible(False)
        ax_histy.spines.top.set_visible(False)


    def get_elbow_dict(
            self,
            method,
            dataset,
            tracks_path,
            cm_size_filter=False,
            cm_len_filter=False,
            min_cm_size=None,
            min_cm_len=None
        ):
        if self.pp_threshold is not None:
            str_ = 'pp_threshold_' + str(self.pp_threshold)
        else:
            str_ = ''
        tracks_output_path = os.path.join(
            self.input_path,
            dataset,
            self.score_type,
            'subset_tracks_and_content',
            self.folder_name
        )
        all_tracks = pd.read_csv(tracks_path, sep='\t', header=None)
        # print(sorted(all_tracks.iloc[:, 2] - all_tracks.iloc[:, 1]))
        if cm_size_filter:
            all_tracks = all_tracks[all_tracks.iloc[:, 9] >= min_cm_size].copy()
        if cm_len_filter:
            all_tracks.loc[:, 'length'] = all_tracks.iloc[:, 2] - all_tracks.iloc[:, 1]
            all_tracks = all_tracks[all_tracks.loc[:, 'length'] >= min_cm_len].copy()
        n_all_CMs = all_tracks.shape[0]
        nCMs_per_threshold_dict = {}
        for threshold in self.similarity_threshold:
            if threshold == 0:
                n_CMs = n_all_CMs
            else:
                tracks_df = self.load_tracks(
                    method,
                    dataset,
                    np.round(threshold, 3),
                    tracks_output_path
                    )
                if cm_size_filter:
                    tracks_df = tracks_df[tracks_df.iloc[:, 9] >= min_cm_size].copy()
                if cm_len_filter:
                    tracks_df.loc[:, 'length'] = tracks_df.iloc[:, 2] - tracks_df.iloc[:, 1]
                    tracks_df = tracks_df[tracks_df.loc[:, 'length'] >= min_cm_len].copy()
                n_CMs = tracks_df.shape[0]
            nCMs_per_threshold_dict[threshold] = ( n_CMs / n_all_CMs ) * 100
        return nCMs_per_threshold_dict

    def elbow_plot(
            self,
            dataset,
            tracks_paths_dict,
            method_color_dict,
            cm_size_filter=False,
            cm_len_filter=False,
            min_cm_size=None,
            min_cm_len=None
        ):
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        for i, method in enumerate(self.methods_list):
            n_modules_per_threshold_dict = self.get_elbow_dict(
                method,
                dataset,
                tracks_paths_dict[method],
                cm_size_filter=cm_size_filter,
                cm_len_filter=cm_len_filter,
                min_cm_size=min_cm_size,
                min_cm_len=min_cm_len
                )
            ax.plot(
                n_modules_per_threshold_dict.keys(),
                n_modules_per_threshold_dict.values(),
                'o-',
                label=method,
                color=method_color_dict[method]
            )

        handles, labels = plt.gca().get_legend_handles_labels()
        ax.legend(handles, labels)
        ax.set_xlabel('Similarity score (ss.)')
        ax.set_ylabel('% chromatin modules across all datasets > ss.')
        ax.set_xlim((0, 1))
        ax.set_ylim((0, 105))
        ax.set_title(dataset)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        if self.save_figure:
            if self.pp_threshold is not None:
                str_ = 'pp_threshold_' + str(self.pp_threshold)
            else:
                str_ = ''
            figure_out_path = os.path.join(
                self.input_path,
                'plots',
                dataset,
                self.folder_name,
                self.score_type + str_,
            )
            if not os.path.exists(figure_out_path):
                os.makedirs(figure_out_path)
            if self.use_cmQTLs:
                plt.savefig(
                    os.path.join(
                        figure_out_path,
                        '_'.join([dataset] + self.methods_list + ['with_cmQTLs_elbow.pdf'])
                        ),
                    dpi=300,
                    bbox_inches='tight'
                    )
            else:
                plt.savefig(
                    os.path.join(
                        figure_out_path,
                        '_'.join([dataset] + self.methods_list + ['all_CMs_elbow.pdf'])
                        ),
                    dpi=300,
                    bbox_inches='tight'
                    )

    def scores_hist_plot(
            self,
            dataset,
            method_A,
            method_B,
            method_color_dict
        ):
        if self.pp_threshold is not None:
            str_ = 'pp_threshold_' + str(self.pp_threshold)
        else:
            str_ = ''
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        A_B_F1_scores = self.load_npy(
            method_A,
            method_B,
            dataset,
            os.path.join(
                self.input_path,
                dataset,
                self.score_type,
                self.score_type + '_scores'
            )
        )
        sim_scores_df = pd.DataFrame(
            [
                [method_A + '_' + cm_ids_A_B.split('_')[0],
                method_B + '_' + cm_ids_A_B.split('_')[1],
                sim_score]
                for cm_ids_A_B, sim_score in A_B_F1_scores.items()
            ],
            columns=['cms_A', 'cms_B', 'score']
        )
        sim_scores_df = sim_scores_df[~sim_scores_df.loc[:, 'cms_A'].str.contains('NAN')]
        sim_scores_df = sim_scores_df[~sim_scores_df.loc[:, 'cms_B'].str.contains('NAN')]
        ax.hist(
            sim_scores_df.loc[:, 'score'],
            bins=50,
            density=True,
            label=method_A + ' vs ' + method_B,
            color='lightgray'
        )
        ax.axvline(
            np.mean(
                sim_scores_df.loc[:, 'score']
            ),
            linestyle='--',
            color='gray',
            linewidth=1
            )
        # ax.axvline(
        #     np.mean(
        #         scores_list[1]
        #     ),
        #     linestyle='--',
        #     color=method_color_dict[methods_[1]],
        #     linewidth=1
        #     )
        ax.set_xlabel('Similarity score (ss.)')
        ax.set_ylabel('Density')
        ax.set_xlim((0, 1))
        # ax.set_ylim((0, 105))
        ax.set_title(dataset)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        if self.save_figure:
            if self.pp_threshold is not None:
                str_ = 'pp_threshold_' + str(self.pp_threshold)
            else:
                str_ = ''
            figure_out_path = os.path.join(
                self.input_path,
                'plots',
                dataset,
                self.folder_name,
                self.score_type + str_,
            )
            if not os.path.exists(figure_out_path):
                os.makedirs(figure_out_path)
            if self.use_cmQTLs:
                plt.savefig(
                    os.path.join(
                        figure_out_path,
                        '_'.join([dataset] + self.methods_list + ['with_cmQTLs_hist.pdf'])
                        ),
                    dpi=300,
                    bbox_inches='tight'
                    )
            else:
                plt.savefig(
                    os.path.join(
                        figure_out_path,
                        '_'.join([dataset] + self.methods_list + ['all_CMs_hist.pdf'])
                        ),
                    dpi=300,
                    bbox_inches='tight'
                    )

    # def enrichr_heatmap(
    #         self,
    #         enrichr_df,
    #         figsize=(5, 5),
    #         colors_dict=None,
    #         title=None,
    #         figure_out_path=None
    #     ):
    #     plt.figure(figsize=(2, 5))
    #     sns.clustermap(
    #         enrichr_df[['CRD', 'VCM', 'PHM']].astype(float),
    #         row_colors=enrichr_df['dataset'].map(colors_dict),
    #         figsize=figsize,
    #         row_cluster=False,
    #         square=True,
    #         yticklabels=True,
    #         cbar_kws=dict(
    #             label='$-log_{10}(~adjusted~pv~)$',
    #             pad=0.01,
    #             shrink=0.3),
    #         mask=(enrichr_df[['CRD', 'VCM', 'PHM']].astype(float)==0.0))

    #     handles = [patches.Patch(facecolor=colors_dict[name]) for name in colors_dict.keys()]
    #     plt.legend(handles,
    #             colors_dict,
    #             title='Cell types',
    #             bbox_to_anchor=(1, 1),
    #             bbox_transform=plt.gcf().transFigure,
    # #                frameon=False,
    #             loc='upper right')
    #     if self.save_figure:
    #         if self.use_cmQTLs:
    #             plt.savefig(
    #                 os.path.join(
    #                     figure_out_path,
    #                     '_'.join([dataset] + self.methods_list  + ['with_cmQTLs', title, 'enrichr_heatmap.pdf'])
    #                     ),
    #                 dpi=300,
    #                 bbox_inches='tight'
    #                 )
    #         else:
    #             plt.savefig(
    #                 os.path.join(
    #                     figure_out_path,
    #                     '_'.join([dataset] + self.methods_list  + ['all_CMs', title, 'enrichr_heatmap.pdf'])
    #                     ),
    #                 dpi=300,
    #                 bbox_inches='tight'
    #                 )

    def enrichr_dotplot(
            self,
            enrichr_df,
            similarity_threshold,
            p_value_col,
            gene_annotation_col,
            n_gene_hits,
            figsize=(5, 5),
            alpha=1,
            set_title=False,
            colors_dict=None,
            title=None,
            figure_name=None
        ):
        similarity_threshold = np.round(similarity_threshold, 3)
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        sp = plt.scatter(
            enrichr_df.loc[:, p_value_col],
            enrichr_df.loc[:, gene_annotation_col],
            alpha=alpha,
            s=5 * np.array(list(enrichr_df.loc[:, n_gene_hits])),
            c=enrichr_df.loc[:, 'dataset'].map(colors_dict)
            )
        ax.set_xlabel('$-log_{10}(~adjusted~pv~)$')
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        if set_title:
            ax.set_title(title)
        ax.set_ylim((-1, len(set(enrichr_df.loc[:, gene_annotation_col]))))
        ax.set_xlim((1, max(enrichr_df.loc[:, p_value_col]) + 0.5))
        legend_1 = ax.legend(
            *sp.legend_elements('sizes', num=3),
            title='Gene set',
            title_fontsize=12,
            loc='upper right',
            bbox_to_anchor=(1.5, 1)
            )
        handles = [patches.Patch(facecolor=colors_dict[name]) for name in colors_dict.keys()]
        legend_2 = ax.legend(
            handles,
            colors_dict,
            title='Cell types',
            title_fontsize=12,
            loc='upper right',
            bbox_to_anchor=(2.1, 1))
        ax.add_artist(legend_1)
        ax.add_artist(legend_2)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        if self.save_figure:
            if self.pp_threshold is not None:
                str_ = '/pp_threshold_' + str(self.pp_threshold)
            else:
                str_ = ''
            figure_out_path = os.path.join(
                self.input_path,
                'plots',
                'enrichr_dotplots',
                self.folder_name,
                self.score_type + str_,
            )
            if not os.path.exists(figure_out_path):
                os.makedirs(figure_out_path)
            if figure_name is not None:
                plt.savefig(
                    os.path.join(
                        figure_out_path,
                        figure_name
                        ),
                    dpi=300,
                    bbox_inches='tight'
                    )
            else:
                if self.use_cmQTLs:
                    plt.savefig(
                        os.path.join(
                            figure_out_path,
                            '_'.join(
                                ['all_datasets'] + self.methods_list +
                                ['only_CMs_with_cmQTL', title, 'enrichr', 'dotplot', str(similarity_threshold) + '.pdf']
                            )
                        ),
                        dpi=300,
                        bbox_inches='tight'
                        )
                else:
                    plt.savefig(
                        os.path.join(
                            figure_out_path,
                            '_'.join(
                                ['all_datasets'] + self.methods_list +
                                ['all_CMs', title, 'enrichr', 'dotplot', str(similarity_threshold) + '.pdf']
                            )
                        ),
                        dpi=300,
                        bbox_inches='tight'
                        )
