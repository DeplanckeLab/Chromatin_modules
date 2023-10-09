import argparse
import logging
import os

from edgelist2VCMs import edgelist2VCMs

parser = argparse.ArgumentParser(
    description='Build VCMs from correlation edgelist ©Olga Pushkarev')
parser.add_argument(
    '-d',
    '--dataset',
    type=str,
    help='Input dataset',
    required=True)
parser.add_argument(
    '-corr_p',
    '--correlation_path',
    type=str,
    help='Path to folder with correlation edgelists',
    required=True)
parser.add_argument(
    '-pv',
    '--pvalue',
    help='p-value threshold as a float or string. To analyze data with multiple p-value thresholds, provide comma-separated input list of p-values',
    required=True)
parser.add_argument(
    '-t',
    '--remove_totem',
    type=int,
    help='Remove totem VCMs from the analysis',
    required=True)
parser.add_argument(
    '-f',
    '--save_files',
    help='Save VCMs tracks and content files',
    type=int,
    required=True)
parser.add_argument(
    '-m',
    '--save_meta',
    help='Save summary statistics for VCMs into dictionary',
    default=False,
    required=False)
parser.add_argument(
    '-o',
    '--output_path',
    type=str,
    help='Output directory for VCM tracks and content files',
    required=True)

args = parser.parse_args()

# Arguments
dataset = args.dataset
path_to_correlation_edgelist =args.correlation_path
pvalue_list = args.pvalue.split(',')
remove_totem = args.remove_totem
save_files = args.save_files
save_meta = args.save_meta
path_to_output_directory = args.output_path

path_to_correlation_edgelist = path_to_correlation_edgelist.replace('\\', '/')
path_to_output_directory = path_to_output_directory.replace('\\', '/')

if not os.path.exists(path_to_output_directory):
    os.makedirs(path_to_output_directory)

logging.basicConfig(
    filename=os.path.join(path_to_output_directory, dataset + '_logfile.txt'),
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    filemode='a'
)
logging.info('########################################################')
logging.info('# Build VCMs from correlation edgelist ©Olga Pushkarev #')
logging.info('########################################################')
logging.info('\nList of arguments:\n')
logging.info('Input:')
logging.info('\t 1. Input dataset: ' + dataset)
logging.info('\t 2. Using the following p-value thresholds: ' + ', '.join(pvalue_list))
logging.info('\t 3. Path to correlation edgelist: ' + path_to_correlation_edgelist)
logging.info('\t 4. Removing totem VCMs: ' + str(remove_totem))
logging.info('\t 5. Saving VCM tracks and content files: ' + str(save_files))
logging.info('\t 6. Saving meta information for VCMs: ' + str(save_meta))
logging.info('\n')
logging.info('\t 1. Output directory for VCM tracks and content files: ' + path_to_output_directory)
logging.info('\n')

edgelist_to_vcms = edgelist2VCMs(
    output_path=path_to_output_directory,
    remove_totem=remove_totem,
    save_files=save_files,
    save_meta=save_meta
)
all_correlations_df = edgelist_to_vcms.load_corr_edgelist(path_to_correlation_edgelist)
if len(pvalue_list) > 1:
    for pv in pvalue_list:
        corrected_pv_df = edgelist_to_vcms.correct_pv(all_correlations_df, float(pv))

        _, _, meta_dict = edgelist_to_vcms.get_VCMs(
            dataset,
            corrected_pv_df,
            str(pv)
        )
        if save_meta:
            edgelist_to_vcms.save_dict(
                '_'.join(
                    'stats_for_VCMs',
                    str(pv),
                    'treshold.npy'
                    ),
                path_to_output_directory,
                meta_dict
            )
else:
    corrected_pv_df = edgelist_to_vcms.correct_pv(all_correlations_df, float(pvalue_list[0]))
    _, _, meta_dict = edgelist_to_vcms.get_VCMs(
        dataset,
        corrected_pv_df,
        str(pvalue_list[0])
    )
    if save_meta:
        edgelist_to_vcms.save_dict(
            '_'.join([
                dataset,
                'stats_for_VCMs',
                str(pvalue_list[0]),
                'treshold.npy'
                ]),
            dataset,
            meta_dict
        )
