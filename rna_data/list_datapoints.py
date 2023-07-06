
import os
import numpy as np
import json
import re
from .datapoint import Datapoint, DatapointFactory
from typing import List, Tuple, Union, Optional
from .parsers import Fasta, DreemOutput

class ListofDatapoints:

    def __init__(self, datapoints):
        self.datapoints = datapoints

    def __call__(self) -> List[Datapoint]:
        return self.datapoints

    @classmethod
    def from_fasta(cls, fasta_file, predict_structure, predict_dms):
        """Create a list of datapoint from a fasta file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        sequences, references = Fasta.parse(fasta_file)
        return cls([DatapointFactory.from_fasta(sequence, reference, predict_structure, predict_dms)
                    for sequence, reference in zip(sequences, references)])

    @classmethod
    def from_ct(cls, ct_files, predict_dms):
        """Create a list of datapoint from a ct file. The dms will be predicted if predict_dms is True."""
        return cls([DatapointFactory.from_ct(ct_file, predict_dms) for ct_file in ct_files])

    @classmethod
    def from_dreem_output(cls, dreem_output_file, predict_structure):
        """Create a list of datapoint from a dreem output file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        return cls([DatapointFactory.from_dreem_output(reference, sequence, mutation_rate, predict_structure)
            for reference, sequence, mutation_rate in DreemOutput.parse(dreem_output_file)])

    @classmethod
    def from_json(cls, json_file):
        """Create a list of datapoint from a json file."""
        with open(json_file) as f:
            return cls([DatapointFactory.from_json_line(line) for line in f if line.strip() not in ['{', '}']])


    def to_base_pairs_npy(self, path):
        """Writes the structure npy file from the list of datapoints.

        Args:
            path (str): path to the npy file

        Returns:
            np.array: structure matrix of pairs of 0-based indices.

        Examples:
            >>> datapoints = ListofDatapoints([Datapoint(reference='reference', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0])])
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
            >>> datapoints = ListofDatapoints([Datapoint(reference='reference', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0])])
            >>> datapoints.to_dms_npy('temp/dms.npy')
            array([[1.0, 2.0, 3.0]], dtype=object)
        """

        arr = np.array([datapoint.dms for datapoint in self.datapoints], dtype=object)
        np.save(path, arr)
        return arr


    def to_sequence_npy(self, path):
        """Writes the sequence npy file from the list of datapoints.

        Args:
            path (str): path to the npy file

        Returns:
            np.array: int-encoded sequence matrix

        Examples:
            >>> datapoints = ListofDatapoints([Datapoint(reference='reference', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0]),\
                                               Datapoint(reference='reference', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0])])
            >>> datapoints.to_sequence_npy('temp/sequence.npy')
            array([[1, 1, 2, 2, 3, 3],
                   [1, 1, 2, 2, 3, 3]], dtype=object)
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
            >>> datapoints = ListofDatapoints([Datapoint(reference='reference', sequence='AACCGG', paired_bases=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0])])
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