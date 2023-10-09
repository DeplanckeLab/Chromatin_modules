import argparse
import logging
import os
import pandas as pd
import scipy.stats as stats

parser = argparse.ArgumentParser(
    description='Build VCMs from correlation edgelist ©Olga Pushkarev')
parser.add_argument(
    '-d',
    '--dataset',
    type=str,
    help='Input dataset',
    required=True)
parser.add_argument(
    '-i',
    '--path_to_input_directory',
    type=str,
    help='Path to folder with theoretical p-values for correlations',
    required=True)
parser.add_argument(
    '-o',
    '--output_path',
    type=str,
    help='Output directory for correlation values and empirical p-values',
    required=True)

args = parser.parse_args()

# Arguments
dataset = args.dataset
path_to_input_directory = args.path_to_input_directory
path_to_output_directory = args.output_path

path_to_input_directory = path_to_input_directory.replace('\\', '/')
path_to_output_directory = path_to_output_directory.replace('\\', '/')

if not os.path.exists(path_to_output_directory):
    os.makedirs(path_to_output_directory)

logging.basicConfig(
    filename=os.path.join(path_to_output_directory, dataset + '_logfile.txt'),
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    filemode='a'
)
logging.info('###########################################################')
logging.info('# Get empirical p-values for correlations ©Olga Pushkarev #')
logging.info('###########################################################')
logging.info('\nList of arguments:\n')
logging.info('Input:')
logging.info('\t 1. Input dataset: ' + dataset)
logging.info('\t 2. Path to correlation edgelist(s): ' + path_to_input_directory)
logging.info('\n')
logging.info('\t 1. Output directory for correlation files: ' + path_to_output_directory)
logging.info('\n')

def empirical_pvalue_for_corr(dataset,
                              input_path,
                              output_path):
    with open(os.path.join(output_path,
                   dataset + '_all_marks_empirical_corr_p_values.txt'), 'a+') as file_with_corr:
        background_distribution = pd.read_csv(input_path,
                                              sep='\t',
                                              header=None)
        background_distribution.columns = ['peak1', 'peak2', 'corr', 'pvalue']
        background_distribution['peak_pair'] = list(zip(background_distribution['peak1'],
                                                        background_distribution['peak2']))
        
        background_mean = background_distribution['corr'].mean()
        background_std_dev = background_distribution['corr'].std()
        dict_with_corr = dict(zip(background_distribution['peak_pair'],
                                  background_distribution['corr']))
        del background_distribution
        for peak_pair, corr in dict_with_corr.items():
            peak1, peak2 = peak_pair
            p_value = 1 - stats.norm.cdf(corr, background_mean, background_std_dev)
            file_with_corr.write(
                '\t'.join(
                    [
                        peak1,
                        peak2,
                        str(corr),
                        str(p_value)
                    ]
                ) + '\n'
            )

input_path = os.path.join(
    path_to_input_directory,
    dataset + '_all_marks_VCM_theoretical_corr_with_p_values.txt'
)

empirical_pvalue_for_corr(dataset, input_path, path_to_output_directory)
