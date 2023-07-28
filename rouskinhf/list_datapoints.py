
import os
import numpy as np
import json
import re
from .datapoint import Datapoint, DatapointFactory
from typing import List, Tuple, Union, Optional
from .parsers import Fasta, DreemOutput
import pandas as pd
from tqdm import tqdm as tqdm_parser

class ListofDatapoints:

    def __init__(self, datapoints=[], verbose=True):
        self.datapoints, self.filtering_report = self.filter_duplicates_and_unvalid(datapoints, verbose)
        
    def __call__(self) -> List[Datapoint]:
        return self.datapoints

    @classmethod
    def from_fasta(cls, fasta_file, predict_structure, predict_dms, tqdm=True, verbose=True):
        """Create a list of datapoint from a fasta file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        sequences, references = Fasta.parse(fasta_file)
        return cls([DatapointFactory.from_fasta(sequence, reference, predict_structure, predict_dms)
                    for sequence, reference in tqdm_parser(zip(sequences, references), total=len(sequences), desc='Parsing fasta file', disable=not tqdm)], verbose=verbose)

    @classmethod
    def from_ct(cls, ct_files, predict_dms, tqdm=True, verbose=True):
        """Create a list of datapoint from a ct file. The dms will be predicted if predict_dms is True."""
        return cls([DatapointFactory.from_ct(ct_file, predict_dms) for ct_file in tqdm_parser(ct_files, total=len(ct_files), desc='Parsing ct files', disable=not tqdm)], verbose=verbose)

    @classmethod
    def from_dreem_output(cls, dreem_output_file, predict_structure, tqdm=True, verbose=True):
        """Create a list of datapoint from a dreem output file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        n_lines = len(list(DreemOutput.parse(dreem_output_file)))
        return cls([DatapointFactory.from_dreem_output(reference, sequence, mutation_rate, predict_structure)
            for reference, sequence, mutation_rate in tqdm_parser(DreemOutput.parse(dreem_output_file), total=n_lines, desc='Parsing dreem output file', disable=not tqdm)], verbose=verbose)

    @classmethod
    def from_json(cls, json_file, tqdm=True, verbose=True):
        """Create a list of datapoint from a json file."""
        with open(json_file) as f:
            return cls([DatapointFactory.from_json_line(line) for line in tqdm_parser(f, total=sum(1 for line in open(json_file) if line.strip() not in ['{', '}']), desc='Parsing json file', disable= not tqdm) if line.strip() not in ['{', '}']], verbose=verbose)


    def to_base_pairs_npy(self, path):
        """Writes the structure npy file from the list of datapoints.

        Args:
            path (str): path to the npy file

        Returns:
            np.array: structure matrix of pairs of 0-based indices.

        Examples:
            >>> datapoints = ListofDatapoints([Datapoint(reference='reference', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0])], verbose=False)
            >>> datapoints.to_base_pairs_npy('temp/base_pairs.npy')
            array([[[1, 2],
                    [3, 4]]], dtype=object)
        """

        arr = np.array([np.array(datapoint.paired_bases) for datapoint in self.datapoints], dtype=object)
        np.save(path, arr)
        return arr


    def to_dms_npy(self, path):
        """Writes the dms npy file from the list of datapoints.

        Args:
            path (str): path to the npy file

        Returns:
            np.array: dms matrix of floats arrays.

        Examples:
            >>> datapoints = ListofDatapoints([Datapoint(reference='reference', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0])], verbose=False)
            >>> datapoints.to_dms_npy('temp/dms.npy')
            array([[1.0, 2.0, 3.0]], dtype=object)
        """
        
        arr = np.array([np.array(datapoint.dms, dtype=np.float32) for datapoint in self.datapoints], dtype=object)
        np.save(path, arr, allow_pickle=True)
        return arr


    def to_sequence_npy(self, path):
        """Writes the sequence npy file from the list of datapoints.

        Args:
            path (str): path to the npy file

        Returns:
            np.array: int-encoded sequence matrix

        Examples:
            >>> from numpy import array
            >>> datapoints = ListofDatapoints([Datapoint(reference='reference', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0]),\
                                               Datapoint(reference='reference', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0])], verbose=False)
            >>> assert not (datapoints.to_sequence_npy('temp/sequence.npy') - array([[1, 1, 2, 2, 3, 3],[1, 1, 2, 2, 3, 3]], dtype=object)).any(), "The sequence matrix is not correct."
        """

        arr = np.array([datapoint.embed_sequence() for datapoint in self.datapoints], dtype=object)
        np.save(path, arr)
        return arr


    def to_reference_npy(self, path):
        """Writes the reference npy file from the list of datapoints.

        Args:
            path (str): path to the npy file

        Returns:
            np.array: reference matrix of int-encoded.

        Examples:
            >>> datapoints = ListofDatapoints([Datapoint(reference='reference', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0])], verbose=False)
            >>> datapoints.to_reference_npy('temp/reference.npy')
            array(['reference'], dtype='<U9')
        """

        arr = np.array([datapoint.reference for datapoint in self.datapoints])
        np.save(path, arr)
        return arr


    def to_json(self, path) -> None:
        """Write a list of datapoints to a json file."""
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            f.write('{\n')
            for idx, datapoint in enumerate(self.datapoints): # write the datapoints one by one to avoid memory issues
                f.write(str(datapoint)+ ',\n' if idx != len(self.datapoints)-1 else str(datapoint)+'\n')
            f.write('}\n')

    def to_pandas(self, datapoints=None) -> pd.DataFrame:
        """Converts the list of datapoints into a pandas dataframe.

        Example:
        >>> datapoints = ListofDatapoints([Datapoint(reference='reference', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0])], verbose=False)
        >>> datapoints.to_pandas()
           reference sequence      paired_bases              dms
        0  reference   AACCGG  ((1, 2), (3, 4))  (1.0, 2.0, 3.0)
        """
        if datapoints is None: datapoints = self.datapoints
        return pd.DataFrame([datapoint.to_flat_dict() for datapoint in datapoints])
    
    def from_pandas(self, df: pd.DataFrame) -> None:
        """Converts a pandas dataframe into a list of datapoints.

        Example:
            >>> df = pd.DataFrame([{'reference': 'reference', 'sequence': 'AACCGG', 'paired_bases': [[1, 2], [3, 4]], 'dms': [1.0, 2.0, 3.0]}])
            >>> ListofDatapoints(verbose=False).from_pandas(df)
            [Datapoint('reference', sequence='AACCGG', paired_bases=((1, 2), (3, 4)), dms=(1.0, 2.0, 3.0))]
        """
        return [Datapoint.from_flat_dict(datapoint_dict) for datapoint_dict in df.to_dict('records')]

    def filter_duplicates_and_unvalid(self, datapoints, verbose=True):
        """Filters out duplicate sequences.
        Only keep the first occurence of a sequence if all the other structures are the same.

        Examples:
            >>> datapoints = ListofDatapoints([ Datapoint(reference='ref1', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[0,0,0,0,0,0]),\
                                Datapoint(reference='ref2', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[0,0,0,0,0,0]),\
                                Datapoint(reference='ref2', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[0,0,0,0,0,0]),\
                                Datapoint(reference='ref4', sequence='AUGGC', paired_bases=[[1, 2]], dms=[0,0,0,0,1]),\
                                Datapoint(reference='ref5', sequence='AUGGC', paired_bases=[[0, 4], [2, 3]], dms=[0,0,0,0,0]),\
                                Datapoint(reference='ref6', sequence='not a regular sequence', paired_bases=[[0, 4]], dms=[0,0,0,0,0]) ])
            Over a total of 6 datapoints, there are:
                - 1 valid datapoints
                - 1 invalid datapoints (ex: sequence with non-regular characters)
                - 1 datapoints with the same reference
                - 1 duplicate sequences with the same structure / dms
                - 2 duplicate sequences with different structure / dms
            >>> datapoints.datapoints
            [Datapoint('ref1', sequence='AACCGG', paired_bases=((1, 2), (3, 4)), dms=(0, 0, 0, 0, 0, 0))]
        """
        
        # Remove None datapoints        
        n_input_datapoints = len(datapoints)
        datapoints = list(filter(lambda x: x is not None, datapoints))
        n_unvalid_datapoints = n_input_datapoints - len(datapoints)

        # Convert to pandas
        df = self.to_pandas(datapoints)
        
        # Remove duplicate or conflicting datapoints and keep track of the number of datapoints removed
        def drop_duplicates(df:pd.DataFrame, **kwargs):
            len_df_before = len(df)
            df.drop_duplicates(**kwargs)
            return len_df_before - len(df)

        
        # Don't allow multiple sequences with the same reference
        n_same_ref_datapoints = drop_duplicates(df, subset=['reference'], inplace=True)

        # Keep only one datapoint per sequence and structure
        if 'paired_bases' in df.columns:
            n_duplicates_datapoints = drop_duplicates(df, subset=['sequence','paired_bases'], inplace=True, ignore_index=True)

        # Keep only one datapoint per sequence and dms
        if 'dms' in df.columns:
            if not 'paired_bases' in df.columns: n_duplicates_datapoints = 0
            n_duplicates_datapoints += drop_duplicates(df, subset=['sequence','dms'], inplace=True, ignore_index=True)
        
        # If there are multiple structures / dms with the same structure, keep none
        n_same_seq_datapoints = drop_duplicates(df, subset=['sequence'], inplace=True, ignore_index=True, keep='first' if not ('paired_bases' in df.columns or 'dms' in df.columns) else False)

        # Convert back to list of datapoints
        datapoints = self.from_pandas(df)

        # Write report
        report = f"""Over a total of {n_input_datapoints} datapoints, there are:
    - {len(datapoints)} valid datapoints
    - {n_unvalid_datapoints} invalid datapoints (ex: sequence with non-regular characters)
    - {n_same_ref_datapoints} datapoints with the same reference"""
        if 'paired_bases' in df.columns or 'dms' in df.columns:
            report += f"""
    - {n_duplicates_datapoints} duplicate sequences with the same structure / dms
    - {n_same_seq_datapoints} duplicate sequences with different structure / dms"""
        else:
            report += f"""
    - {n_same_seq_datapoints} duplicate sequences"""
    
        if verbose: print(report)
            
        return datapoints, report
    
            
def base_pairs_into_matrix(base_pairs):
    """Turns a list of base pairs into a matrix of 1s and 0s.

    Args:
        base_pairs (list): list of base pairs, where each base pair is a tuple of 0-based two ints.

    Returns:
        matrix (np.array): matrix of 1s and 0s.

    Example:
        >>> from numpy import array, int64
        >>> base_pairs = [(0, 1), (2, 3), (4, 5)]
        >>> assert (base_pairs_into_matrix(base_pairs) ==  array([[0., 1., 0., 0., 0., 0.], \
                                                                [1., 0., 0., 0., 0., 0.], \
                                                                [0., 0., 0., 1., 0., 0.], \
                                                                [0., 0., 1., 0., 0., 0.], \
                                                                [0., 0., 0., 0., 0., 1.], \
                                                                [0., 0., 0., 0., 1., 0.]], dtype=int64)).all()

    """
    max_index = max([max(base_pair) for base_pair in base_pairs])
    arr = np.zeros((max_index + 1, max_index + 1))
    for base_pair in base_pairs:
        arr[base_pair[0], base_pair[1]] = 1
        arr[base_pair[1], base_pair[0]] = 1
    return arr