import datetime
import itertools
import networkx as nx
import numpy as np
import operator
import os
import pandas as pd
import rpy2.robjects as ro
import signal
import sys
sys.path.append('/usr/local/lib/python3.9/site-packages')

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException

# Change the behavior of SIGALRM
signal.signal(signal.SIGALRM, timeout_handler)

def solve_cyclic_graph(G, timeout, j=0):
    G_copy = G.copy()
    G_edges_dict = {
        (edge[0], edge[1]): edge[2]['weight']
        for edge in G.edges(data=True)
    }
    G_edges_dict_sorted = {
        key: G_edges_dict[key]
        for key in sorted(
            G_edges_dict,
            key=G_edges_dict.get,
            reverse=False
        )
    }
    G.remove_edge(
        *min(
            G_edges_dict.items(),
            key=operator.itemgetter(1)
        )[0]
    )
    if nx.is_weakly_connected(G) and nx.is_directed_acyclic_graph(G):
        return G
    elif nx.is_weakly_connected(G) and not nx.is_directed_acyclic_graph(G):
        if datetime.datetime.now() > timeout:
            return None
        return solve_cyclic_graph(G, timeout, j=0)
    else:
        j += 1
        G_copy.remove_edge(*list(G_edges_dict_sorted)[j])
        if nx.is_weakly_connected(G_copy) and nx.is_directed_acyclic_graph(G_copy):
            return G_copy
        elif nx.is_weakly_connected(G_copy) and not nx.is_directed_acyclic_graph(G_copy):
            if datetime.datetime.now() > timeout:
                return None
            return solve_cyclic_graph(G_copy, timeout, j)
        else:
            if datetime.datetime.now() > timeout:
                return None
            return solve_cyclic_graph(G_copy, timeout, j)


def find_the_heaviest_hamiltonian_path(G):
    G_nodes = list(G.nodes())
    G_edges_dict = {(edge[0], edge[1]): edge[2]['weight']
                    for edge in G.edges(data=True)}
    source_dest_combinations = list(itertools.combinations(G_nodes, 2))
    heaviest_path_dict = {'path': None,
                          'path_weight': 0}
    for source, dest in source_dest_combinations:
        for path in nx.all_simple_paths(G, source, dest):
            w = 0
            if len(path) == len(G_nodes):
                for e in zip(path[:-1], path[1:]):
                    w += G_edges_dict[e]
                if w > heaviest_path_dict['path_weight']:
                    heaviest_path_dict['path'] = path
                    heaviest_path_dict['path_weight'] = w
    if heaviest_path_dict['path_weight'] == 0:
        return None
    else:
        heaviest_path = heaviest_path_dict['path']
        g = nx.Graph(((u, v, e) for u, v, e in G.edges_iter(data=True)
                      if (u, v) in zip(heaviest_path[:-1], heaviest_path[1:])))
        return g


def split_dag_into_subpaths(G):
    roots = []
    leaves = []
    for node in G.nodes():
        if G.in_degree(node) == 0:  # it's a root
            roots.append(node)
        elif G.out_degree(node) == 0:  # it's a leaf
            leaves.append(node)
    output_subpaths = []
    for root in roots:
        for leaf in leaves:
            # Start the timer (15 seconds)
            signal.alarm(15)
            try:
                lst = list(nx.all_simple_paths(G, root, leaf))
                if (lst != []) and (len(lst) <= 1000):
                    output_subpaths.append(lst)
            except TimeoutException:
                # continue the for loop if function A takes more than 15 seconds
                continue
    return output_subpaths


def merge_intersecting_intervals(list_of_pairs):
    '''For the given list of pairs the function finds and merges intersecting intervals.
       The output is the list of intervals of length <= length(input list)'''
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

def build_dags(input_path,
               chromosome,
               threshold,
               dataset,
               save_files=False):
    peak_coordinates = pd.read_csv(
        os.path.join(
            input_path,
            chromosome,
            '_'.join([chromosome, 'peak_coordinates.bed.gz'])
        ),
        sep='\t',
        header=None
    )
    if peak_coordinates.shape[1] == 3:
        peak_coordinates.columns = ['chr', 'start', 'end']
#         peak_coordinates['mark'] = 'ATAC'
    else:
        peak_coordinates.columns = ['chr', 'start', 'end', 'mark']
    peak_coordinates['peak_id'] = np.arange(1, peak_coordinates.shape[0] + 1)
    peak_coordinates['peak_id_str'] = \
        peak_coordinates['mark'] + ':' + \
            'chr' + peak_coordinates['chr'].astype(str) + ':' + \
                peak_coordinates['start'].astype(str) + ':' + \
                    peak_coordinates['end'].astype(str)
    pp = pd.read_csv(
        os.path.join(input_path, chromosome, 'phm_output', 'pp.gz'),
        sep='\t',
        header=None
    )
    pp.columns = ['peak_1', 'peak_2', 'H_0', 'H_11', 'H_12', 'linkage',
                  'pleiotropy', 'causality_1', 'causality_2']
    pp = pp.set_index(['peak_1', 'peak_2'])
    pp_theshold = pp[pp > np.log(threshold)]
    pp_theshold = pp_theshold[['causality_1', 'causality_2']]
    pp_filtered = pp_theshold.dropna(how='all')
    pp_filtered.reset_index(inplace=True)
    if pp_filtered.shape[0] == 0:
        return
    if save_files:
        pp_filtered.to_csv(
            os.path.join(
                input_path,
                chromosome,
                'phm_output',
                '_'.join(['pp_filtered_causality', str(threshold), 'threshold.gz'])
            ),
            sep='\t',
            header=False,
            index=False
        )
    pp_filtered.columns = ['from', 'to', 'causality_1', 'causality_2']
    pp_filtered = pp_filtered.fillna(10)
    pp_filtered['causality_1'] = np.exp(pp_filtered['causality_1'])
    pp_filtered['causality_2'] = np.exp(pp_filtered['causality_2'])
    pp_filtered['causality_1'] = pp_filtered['causality_1'].replace(np.exp(10), 0)
    pp_filtered['causality_2'] = pp_filtered['causality_2'].replace(np.exp(10), 0)
    df_1_ = pp_filtered[['from', 'to', 'causality_1']].copy()
    df_2_ = pp_filtered[['to', 'from', 'causality_2']].copy()
    df_1_.columns = ['from', 'to', 'weight']
    df_2_.columns = ['from', 'to', 'weight']
    edges_causality = pd.concat([df_1_, df_2_], axis=0)
    causality_df = edges_causality[edges_causality['weight'] != 0].copy()
    
    g = nx.from_pandas_edgelist(
        causality_df,
        'from',
        'to',
        'weight',
        create_using=nx.MultiDiGraph()
    )
    list_of_subgraphs = list([g.subgraph(c).copy() for c in nx.weakly_connected_components(g)])
    list_with_dags = []
    list_with_tracks = []
    list_with_contents = []
    for index, subgraph in enumerate(list_of_subgraphs):
        time_out = datetime.datetime.now() + datetime.timedelta(minutes=2)
        if not nx.is_directed_acyclic_graph(subgraph):
            try:
                dag = solve_cyclic_graph(subgraph, time_out)
                subpaths = split_dag_into_subpaths(dag)
                list_with_dags.extend(subpaths)
            except:
                try:
                    edges_to_remove = [
                        (edge[0], edge[1])
                        for edge in subgraph.edges(data=True)
                        if edge[2]['weight'] < threshold + 0.1
                    ]
                    subgraph.remove_edges_from(edges_to_remove)
                    dag = solve_cyclic_graph(subgraph, time_out)
                    subpaths = split_dag_into_subpaths(dag)
                    list_with_dags.extend(subpaths)
                except:
                    continue
        else:
            subpaths = split_dag_into_subpaths(subgraph)
            list_with_dags.extend(subpaths)
        peak_df = peak_coordinates[peak_coordinates['peak_id'].isin(subgraph.nodes())]
        start_end = sorted(list(zip(peak_df['start'], peak_df['end'])), key=lambda x: x[0])
        start = np.array([pair[0] for pair in start_end])
        peak_starts = list(start - min(start_end, key=lambda x: x[0])[0])

        list_with_tracks.append(
            {
                'chr': peak_df['chr'].iloc[0],
                'start': min(peak_df['start']),
                'end': max(peak_df['end']),
                'cm_id': 'phm' + str(index + 1) + '_' + chromosome,
                'number': 1000,
                'strain': '+',
                'start_duplicate': min(peak_df['start']),
                'end_duplicate': max(peak_df['end']),
                'numbers': '0,0,0',
                'cm_size': peak_df.shape[0],
                'peak_length': ','.join(map(str, [j - i for i, j in start_end])),
                'peak_starts': ','.join(map(str, peak_starts)),
                'totem_cm': 1 if len(merge_intersecting_intervals(start_end)) == 1 else 0
            }
        )

        list_with_contents.append(
            {
                'cm_id': 'phm' + str(index + 1) + '_' + chromosome,
                'cm_size': peak_df.shape[0],
                'peaks': ','.join(list(peak_df['peak_id_str']))
            }
        )
    if save_files:

        phm_tracks_bed = pd.DataFrame(list_with_tracks)
        phm_tracks_bed.to_csv(
            os.path.join(
                input_path,
                chromosome,
                '_'.join([dataset, str(threshold), 'phm.tracks.bed'])
            ),
            sep='\t',
            index=False,
            header=False
        )
        del phm_tracks_bed
        phm_content = pd.DataFrame(list_with_contents)
        phm_content.to_csv(
            os.path.join(
                input_path,
                chromosome,
                '_'.join([dataset, str(threshold), 'phm.content.txt'])
            ),
            sep='\t',
            header=False,
            index=False
        )
        del phm_content

    list_with_dags_single_lst = [sublist for lst in list_with_dags for sublist in lst]
    dag_df = pd.DataFrame((el for el in itertools.zip_longest(*list_with_dags_single_lst))).T
    dag_df = dag_df.fillna(0).astype(int)
    if save_files:
        dag_df.to_csv(
            os.path.join(
                input_path,
                chromosome,
                '_'.join(['DAGs_all_chr', str(threshold), 'threshold.csv'])
            ),
            sep='\t',
            header=False,
            index=False
        )

# G = nx.DiGraph()
# G.add_nodes_from([1, 2, 3])
# G.add_edges_from([(1, 2), (1, 3), (3, 4), (3, 5)])
# print(nx.is_directed_acyclic_graph(G))
# print(list(nx.topological_sort(G)))

# roots = []
# leaves = []
# for node in G.nodes():
#     if G.in_degree(node) == 0 : # it's a root
#         roots.append(node)
#     elif G.out_degree(node) == 0: # it's a leaf
#         leaves.append(node)
# [path for root in roots for leaf in leaves for path in nx.all_simple_paths(G, root, leaf)]
