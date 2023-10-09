import argparse
import logging
import numpy as np
import os
import tqdm.auto as tqdm

from BuildDAGs_functions import build_dags

parser = argparse.ArgumentParser(
    description='Build DAGs from PHM output ©Olga Pushkarev')
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
    help='Path to folder with count matrices',
    default=0,
    required=True)
parser.add_argument(
    '-t',
    '--pp_threshold',
    type=str,
    help='Comma separated posterior probability thresholds for causal peak interaction hypothesis, e.g., 0.5,0.6,0.7 or 0.5',
    default='0.5',
    required=False)
parser.add_argument(
    '-c',
    '--chromosomes',
    type=str,
    help='Chromosomes to process',
    required=True)

args = parser.parse_args()

# Arguments
dataset = args.dataset
path_to_input_directory = args.path_to_input_directory
pp_threshold = float(args.pp_threshold)
chromosomes_str = args.chromosomes
path_to_output_directory = path_to_input_directory

if ',' in chromosomes_str:
    if not 'chr' in chromosomes_str.split(',')[0]:
        chromosomes = [
            'chr' + str(chr_id)
            for chr_id in chromosomes_str.split(',')
        ]
    else:
        chromosomes = chromosomes_str.split(',')
elif '-' in chromosomes_str:
    chromosomes_lst = [el.replace('chr', '') for el in chromosomes_str.split('-')]
    chromosomes = [
        'chr' + str(chr_int)
        for chr_int in np.arange(int(chromosomes_lst[0]), int(chromosomes_lst[1]) + 1)]
else:
    if not 'chr' in chromosomes_str:
        chromosomes = ['chr' + chromosomes_str]
    else:
        chromosomes = [chromosomes_str]

logging.basicConfig(
    filename=os.path.join(path_to_output_directory, dataset + '_logfile.txt'),
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    filemode='a'
)
logging.info('##############################################')
logging.info('# Build DAGs from PHM output ©Olga Pushkarev #')
logging.info('##############################################')
logging.info('\nList of arguments:\n')
logging.info('Input:')
logging.info('\t 1. Input dataset: ' + dataset)
path_to_input_directory = path_to_input_directory.replace('\\', '/')
logging.info('\t 2. Path to input files): ' + path_to_input_directory)
logging.info('\n')
logging.info('Output:')
logging.info('\t 1. Path to output files: ' + path_to_input_directory)

for chromosome in tqdm.tqdm(chromosomes):
    build_dags(
        path_to_output_directory,
        chromosome,
        pp_threshold,
        dataset,
        save_files=True
    )

logging.info('Next step: merge DAGs, track and content files across all chromosomes with 2.get_VCMs_from_DAGs.py')
