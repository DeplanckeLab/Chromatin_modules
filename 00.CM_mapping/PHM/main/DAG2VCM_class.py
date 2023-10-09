import os
import pandas as pd
import sys
sys.path.append('/usr/local/lib/python3.9/site-packages')

class DAG2VCM:
    def __init__(self,
                 core_path,
                 chromosomes,
                 threshold,
                 dataset):
        self.core_path = core_path
        self.chromosomes = chromosomes
        self.threshold = threshold
        self.dataset = dataset

    def merge_tracks(self):
        '''Merge PHM tracks for different chromosomes'''
        phm_tracks_to_merge = []
        for chromosome in self.chromosomes:
            try:
                phm_tracks_to_merge.append(
                    pd.read_csv(
                        os.path.join(
                            self.core_path,
                            chromosome,
                            '_'.join([self.dataset, str(self.threshold), 'phm.tracks.bed'])
                            ),
                            sep='\t',
                            header=None
                        ).reset_index(drop=True)
                    )
            except:
                continue
        if not phm_tracks_to_merge:
            print('No causal interactions found')
            return pd.DataFrame()
        else:
            return pd.concat(phm_tracks_to_merge, axis=0)

    def merge_contents(self):
        '''Merge PHM contents for different chromosomes'''
        phm_content_to_merge = []
        for chromosome in self.chromosomes:
            try:
                phm_content_to_merge.append(
                    pd.read_csv(
                        os.path.join(
                            self.core_path,
                            chromosome,
                            '_'.join([self.dataset, str(self.threshold), 'phm.content.txt'])
                            ),
                            sep='\t',
                            header=None
                        ).reset_index(drop=True)
                    )
            except:
                continue
        if not phm_content_to_merge:
            print('No causal interactions found')
            return pd.DataFrame()
        else:
            return pd.concat(phm_content_to_merge, axis=0)

    def merge_dags(self):
        '''Merge DAGs from different chromosomes'''
        dags_to_merge = []
        for chromosome in self.chromosomes:
            try:
                dags_chr = pd.read_csv(
                    os.path.join(
                        self.core_path,
                        chromosome,
                        '_'.join(['DAGs_all_chr', str(self.threshold), 'threshold.csv'])
                        ),
                        sep='\t',
                        header=None
                    ).reset_index(drop=True)
                dags_chr = dags_chr.astype(str) + '_' + chromosome
                dags_chr = dags_chr.replace({'0_' + chromosome: '0'})
                dags_to_merge.append(dags_chr)
            except:
                continue
        if not dags_to_merge:
            print('No causal interactions found')
            return pd.DataFrame()
        else:
            return pd.concat(dags_to_merge, axis=0).fillna(0)
