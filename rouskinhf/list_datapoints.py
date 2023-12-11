
import os
import numpy as np
from .datapoint import Datapoint, DatapointFactory
from typing import List
from .parsers import Fasta, DreemOutput
from .util import UKN

import pandas as pd
from sklearn.metrics import roc_auc_score
from tqdm import tqdm as tqdm_parser
import json

class ListofDatapoints:

    def __init__(self, datapoints=[], verbose=True):
        self.datapoints = datapoints
        
    def __call__(self) -> List[Datapoint]:
        return self.datapoints

    @classmethod
    def from_fasta(cls, fasta_file, predict_structure, tqdm=True, verbose=True)->'ListofDatapoints':
        """Create a list of datapoint from a fasta file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        sequences, references = Fasta.parse(fasta_file)
        return cls([DatapointFactory.from_fasta(sequence, reference, predict_structure)
                    for sequence, reference in tqdm_parser(zip(sequences, references), total=len(sequences), desc='Parsing fasta file', disable=not tqdm)], verbose=verbose)
    
    @classmethod
    def from_bpseq(cls, bpseq_files, tqdm=True, verbose=True):
        """Create a list of datapoint from a bpseq file. The dms will be predicted if predict_dms is True."""
        return cls([DatapointFactory.from_bpseq(ct_file) for ct_file in tqdm_parser(bpseq_files, total=len(bpseq_files), desc='Parsing bpseq files', disable=not tqdm)], verbose=verbose)

    @classmethod
    def from_ct(cls, ct_folder, tqdm=True, verbose=True):
        """Create a list of datapoint from a ct file. The dms will be predicted if predict_dms is True."""
        ct_files = [os.path.join(ct_folder, ct_file) for ct_file in os.listdir(ct_folder) if ct_file.endswith('.ct')]
        return cls([DatapointFactory.from_ct(ct_file) for ct_file in tqdm_parser(ct_files, total=len(ct_files), desc='Parsing ct files', disable=not tqdm)], verbose=verbose)

    @classmethod
    def from_dreem_output(cls, dreem_output_file, predict_structure, tqdm=True, verbose=True):
        """Create a list of datapoint from a dreem output file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        n_lines = len(list(DreemOutput.parse(dreem_output_file)))
        return cls([DatapointFactory.from_dreem_output(reference, sequence, mutation_rate, predict_structure)
            for reference, sequence, mutation_rate in tqdm_parser(DreemOutput.parse(dreem_output_file), total=n_lines, desc='Parsing dreem output file', disable=not tqdm)], verbose=verbose)

    @classmethod
    def from_json(cls, json_file, predict_structure=False, tqdm=True, verbose=True):
        """Create a list of datapoint from a json file."""
        data = json.load(open(json_file))
        return cls([DatapointFactory.from_json_line(reference, line, predict_structure) for reference, line in tqdm_parser(data.items(), total=len(data), desc='Parsing json file', disable=not tqdm)], verbose=verbose)

    def to_dict(self) -> dict:
        """Converts the list of datapoints into a dictionary."""
        return {datapoint.reference: datapoint.to_dict()[1] for datapoint in self.datapoints}

    def to_json(self, path) -> None:
        """Write a list of datapoints to a json file."""
        if os.path.dirname(path) != '':
            os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            f.write('{\n')
            for idx, datapoint in enumerate(self.datapoints): # write the datapoints one by one to avoid memory issues
                f.write(str(datapoint)+ ',\n' if idx != len(self.datapoints)-1 else str(datapoint)+'\n')
            f.write('}\n')

    def to_pandas(self, datapoints=None) -> pd.DataFrame:
        """Converts the list of datapoints into a pandas dataframe.

        Example:
        >>> datapoints = ListofDatapoints([Datapoint(reference='reference', sequence='AACCGG', structure=[[1, 2], [3, 4]], dms=[1,0,0,0,0,1])], verbose=False)
        >>> temp = datapoints.to_pandas()
        >>> temp['reference'].iloc[0]
        'reference'
        >>> temp['sequence'].iloc[0]
        'AACCGG'
        >>> temp['structure'].iloc[0]
        ((1, 2), (3, 4))
        >>> temp['dms'].iloc[0]
        (1, 0, 0, 0, 0, 1)
        """
        if datapoints is None: datapoints = self.datapoints
        for dp in datapoints:
            dp.convert_structure_to_tuple()
        return pd.DataFrame([datapoint.to_flat_dict() for datapoint in datapoints])
    
    def from_pandas(self, df: pd.DataFrame) -> None:
        """Converts a pandas dataframe into a list of datapoints.

        Example:
            >>> df = pd.DataFrame([{'reference': 'reference', 'sequence': 'AACCGG', 'structure': [[1, 2], [3, 4]], 'dms': [1.0, 2.0, 3.0]}])
            >>> ListofDatapoints(verbose=False).from_pandas(df)
            [Datapoint('reference', sequence='AACCGG', structure=((1, 2), (3, 4)), dms=(1.0, 2.0, 3.0))]
        """
        df.rename(columns={'structure': 'paired_bases'}, inplace=True)
        datapoints = [Datapoint.from_flat_dict(datapoint_dict) for datapoint_dict in df.to_dict('records')]
        for dp in datapoints:
            dp.convert_structure_to_list()
        return datapoints
    
    def drop_none_dp(self):
        n_input_datapoints = len(self.datapoints)
        self.datapoints = list(filter(lambda x: x is not None, self.datapoints))
        n_unvalid_datapoints = n_input_datapoints - len(self.datapoints)
        return self.datapoints, n_unvalid_datapoints

    def filter(self, min_AUROC:int=0.8):
        """Filters out duplicate sequences.
        Only keep the first occurence of a sequence if all the other structures are the same.

        Examples:
            >>> datapoints = ListofDatapoints([ Datapoint(reference='ref1', sequence='AACCGG', structure=[[1, 2], [3, 4]], dms=[1,0,0,0,0,1]),\
                                Datapoint(reference='ref2', sequence='AACCGG', structure=[[1, 2], [3, 4]], dms=[1,0,0,0,0,1]),\
                                Datapoint(reference='ref2', sequence='AACCGG', structure=[[1, 2], [3, 4]], dms=[1,0,0,0,0,1]),\
                                Datapoint(reference='ref4', sequence='AUGGC', structure=[[1, 2]], dms=[0,0,0,0,1]),\
                                Datapoint(reference='ref5', sequence='AUGGC', structure=[[0, 4], [2, 3]], dms=[0,0,0,0,0]),\
                                Datapoint(reference='ref6', sequence='not a regular sequence', structure=[[0, 4]], dms=[0,0,0,0,0]) ])
            Over a total of 6 datapoints, there are:
                - 1 valid datapoints
                - 1 invalid datapoints (ex: sequence with non-regular characters)
                - 1 datapoints with the same reference
                - 1 duplicate sequences with the same structure / dms
                - 2 duplicate sequences with different structure / dms
                - 0 datapoints removed because of low AUROC (<0.8)
            >>> datapoints.datapoints
            [Datapoint('ref1', sequence='AACCGG', structure=((1, 2), (3, 4)), dms=(1, 0, 0, 0, 0, 1))]
        """
        
        # Remove None datapoints   
        n_input_datapoints = len(self.datapoints)   
        datapoints, n_unvalid_datapoints = self.drop_none_dp()
        
        # Remove bad structures
        bad_structures_idx = []
        for idx, datapoint in enumerate(datapoints):
            if not datapoint._assert_paired_bases():
                bad_structures_idx.append(idx)
        datapoints = [datapoint for idx, datapoint in enumerate(datapoints) if idx not in bad_structures_idx]
        n_bad_structures_datapoints = len(bad_structures_idx)

        # Convert to pandas
        df = self.to_pandas(datapoints)

        # Remove duplicate or conflicting datapoints and keep track of the number of datapoints removed
        def drop_duplicates(df:pd.DataFrame, **kwargs):
            len_df_before = len(df)
            df.drop_duplicates(**kwargs)
            return len_df_before - len(df)
        
        # If multiple sequences with the same reference, rename the reference
        refs = dict()
        n_same_ref_datapoints = 0
        for idx, row in df.iterrows():
            if row['reference'] in refs:
                df.at[idx, 'reference'] = f"{row['reference']}_{refs[row['reference']]}"
                n_same_ref_datapoints += 1
            refs[row['reference']] = 1

        # Keep only one datapoint per sequence and structure
        if 'structure' in df.columns:
            n_duplicates_datapoints = drop_duplicates(df, subset=['sequence','structure'], inplace=True, ignore_index=True, keep='first')

        # Keep only one datapoint per sequence and dms
        if 'dms' in df.columns:
            if not 'structure' in df.columns: n_duplicates_datapoints = 0
            n_duplicates_datapoints += drop_duplicates(df, subset=['sequence','dms'], inplace=True, ignore_index=True, keep='first')
        
        # If there are multiple structures / dms with the same sequence, keep none
        n_same_seq_datapoints = drop_duplicates(df, subset=['sequence'], inplace=True, ignore_index=True, keep='first')

        ## Filter out references with low AUROC 
        if 'structure' in df.columns and 'dms' in df.columns:
            def calculate_auroc(row):
                dms = np.array(row['dms'])
                isUnpaired = np.ones_like(dms)
                isUnpaired[np.array(row['structure']).flatten()] = 0
                if set(isUnpaired[dms!=UKN]) != set([0,1]):
                    return 0
                return roc_auc_score(isUnpaired[dms!=UKN], dms[dms!=UKN])

            # Create a boolean mask for rows with auroc score greater than or equal to a threshold
            mask_high_AUROC = df.apply(lambda row: calculate_auroc(row) >= min_AUROC, axis=1)
            df = df[mask_high_AUROC]

        # Convert back to list of datapoints
        self.datapoints = self.from_pandas(df)

        # Write report
        report = f"""Over a total of {n_input_datapoints} datapoints, there are:
    OUTPUT
    - {len(self.datapoints)} valid datapoints
    MODIFIED
    - {n_same_ref_datapoints} multiple sequences with the same reference (renamed reference)
    FILTERED OUT
    - {n_unvalid_datapoints} invalid datapoints (ex: sequence with non-regular characters)
    - {n_bad_structures_datapoints} datapoints with bad structures"""
        if 'structure' in df.columns or 'dms' in df.columns:
            report += f"""
    - {n_duplicates_datapoints} duplicate sequences with the same structure / dms
    - {n_same_seq_datapoints} duplicate sequences with different structure / dms"""
        else:
            report += f"""
    - {n_same_seq_datapoints} duplicate sequences"""
        if 'structure' in df.columns and 'dms' in df.columns:
            report += f"""
    - {np.sum(~mask_high_AUROC)} datapoints removed because of low AUROC (<{min_AUROC})"""

        return report
    