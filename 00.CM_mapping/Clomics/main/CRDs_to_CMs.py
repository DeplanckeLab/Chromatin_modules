import argparse
import logging
import numpy as np
import os
import pandas as pd
import tqdm
import sys

parser = argparse.ArgumentParser(
    description='Convert CRDs to CMs ©Olga Pushkarev')
parser.add_argument(
    '-d',
    '--dataset',
    type=str,
    help='Input dataset',
    required=True)
parser.add_argument(
    '-i',
    '--path_to_input_file',
    type=str,
    help='Path to folder with correlation edgelists',
    required=True)
parser.add_argument(
    '-o',
    '--output_path',
    type=str,
    help='Output directory for CM tracks and content files',
    required=True)

args = parser.parse_args()

# Arguments
dataset = args.dataset
path_to_input_file =args.path_to_input_file
path_to_output_directory = args.output_path

path_to_output_directory = path_to_output_directory.replace('\\', '/')

logging.basicConfig(
    filename=os.path.join(path_to_output_directory, dataset + '_logfile.txt'),
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    filemode='a'
)
logging.info('########################################')
logging.info('# Convert CRDs to CMs ©Olga Pushkarev #')
logging.info('########################################')
logging.info('\nList of arguments:\n')
logging.info('Input:')
logging.info('\t 1. Input dataset: ' + dataset)
logging.info('\t 2. Input path to Clomics .bed file: ' + path_to_input_file)
logging.info('\n')
logging.info('\t 1. Output directory for CRDs in CM-like format, i.e. tracks and content files: ' + path_to_output_directory)
logging.info('\n')

# Functions
def find_children(lst):
    '''Function that finds children of each node in the input list.
       If the node is not internal (is not epigenetic mark, i.e. H3K4me3, H3K27ac, etc.)
       the functions does not search for its kids and continues to the other nodes.
       If no internal nodes left, the function outputs the list of all children of the nodes in the input list'''
    children_list = []
    for node in lst:
            row = crd_tree_df[crd_tree_df['feature_id'] == str(node)]
            if (str(row['first_child'].iloc[0]) == 'nan') or (str(row['second_child'].iloc[0]) == 'nan'):
                children_list.append(str(row['UID'].iloc[0]))
            else:
                children_list.append(str(row['first_child'].iloc[0]))
                children_list.append(str(row['second_child'].iloc[0]))
        else:
            children_list.append(node)
    return children_list


def find_leafs(list_of_parents):
    '''Recursive function to find marks (leafs) that belong to one CRD'''
    if not any(child.startswith('chr') for child in list_of_parents):
        return list_of_parents
    else:
        children = find_children(list_of_parents)
        return find_leafs(children)


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

# Import data from Clomics output
crd_tree_df = pd.read_csv(path_to_input_file, sep='\t', header=None)
crd_tree_df.columns = ['chr', 'start', 'end', 'UID', 'feature_id', 'first_child', 'second_child', 'bool_CRD']
crd_tree_df = crd_tree_df.astype({'first_child': str, 'second_child': str})

peak_id_to_check = crd_tree_df[crd_tree_df.loc[:, 'first_child'] == 'nan'].loc[0, 'UID'].copy()

if peak_id_to_check.startswith('chr'):
    sys.exit('Error! Rename the peaks to avoid confusion between the names of inner nodes of a tree and its leafs.')
else:
    # Select nodes that were marked as CRDs
    crd_nodes = crd_tree_df[crd_tree_df.loc[:, 'bool_CRD'] == 1]
    crd_nodes.index = np.arange(crd_nodes.shape[0])

    # Create dfs to save CRDs in CM-like format
    bed_df = pd.DataFrame(columns=['chr', 'start', 'end', 'CRD_id', 'number', 'strain',
                                   'start_duplicate', 'end_duplicate', 'numbers', 'CRD_size',
                                   'peak_length', 'peak_starts', 'totem_VCM'],
                          index=np.arange(crd_nodes.shape[0]))
    crd_content = pd.DataFrame(columns=['CRD_id', 'CRD_size', 'peaks'],
                               index=np.arange(crd_nodes.shape[0]))
    for index, crd in tqdm.tqdm(crd_nodes.iterrows()):
        parent1, parent2 = crd['first_child'], crd['second_child']
        list_with_marks = find_leafs([parent1, parent2])  # Find marks (leafs of the tree) that belong to one CRD

        # Fill the fields of the .content file
        crd_content['CRD_id'].loc[index] = 'crd' + str(index)
        crd_content['CRD_size'].loc[index] = int(len(list_with_marks))
        crd_content['peaks'].loc[index] = ','.join(list_with_marks)

        # For marks in CRD prepare list of peak coordinates (start, end)
        start_end = sorted([(int(crd_tree_df[crd_tree_df['UID'] == mark]['start']),
                             int(crd_tree_df[crd_tree_df['UID'] == mark]['end'])) for mark in list_with_marks],
                           key = lambda x: x[0])
        start = np.array([pair[0] for pair in start_end])
        peak_starts = list(start - min(start_end, key = lambda x: x[0])[0])

        # Fill the fields of the .bed file
        bed_df['chr'].loc[index] = crd['chr']
        bed_df['start'].loc[index] = min(start_end, key = lambda x: x[0])[0]
        bed_df['end'].loc[index] = max(start_end, key = lambda x: x[1])[1]
        bed_df['CRD_id'].loc[index] = 'crd' + str(index)
        bed_df['number'].loc[index] = int(1000)
        bed_df['strain'].loc[index] = '+'
        bed_df['start_duplicate'].loc[index] = min(start_end, key = lambda x: x[0])[0]
        bed_df['end_duplicate'].loc[index] = max(start_end, key = lambda x: x[1])[1]
        bed_df['numbers'].loc[index] = '0,0,0'
        bed_df['CRD_size'].loc[index] = int(len(list_with_marks))
        bed_df['peak_length'].loc[index] = ','.join(map(str, [j - i for i, j in start_end]))
        bed_df['peak_starts'].loc[index] = ','.join(map(str, peak_starts))
        bed_df['totem_CM'].loc[index] = 1 if len(merge_intersecting_intervals(start_end)) == 1 else 0

    # Save bed and content files
    bed_df.to_csv(os.path.join(path_to_output_directory, dataset + '_crd.tracks.bed'),
                  sep='\t',
                  index=False,
                  header=False)
    crd_content.to_csv(os.path.join(path_to_output_directory, dataset + '_crd.content.txt'),
                       sep='\t',
                       header=False,
                       index=False)
