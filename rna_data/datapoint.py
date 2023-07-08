
from typing import Any, List
from .parsers import *
from .util import add_braces_if_no_braces, dot2int, int2dot, seq2int
import numpy as np

class Datapoint:

    def __init__(self, sequence, reference, structure=None, dms=None, paired_bases=None):
        for attr in [sequence, reference]:
            assert isinstance(attr, str), f"Expected {attr} to be a string, got {type(attr)} instead."
            assert len(attr) > 0, f"Expected {attr} to be non-empty."
        self.reference = reference
        self.sequence = sequence
        self.structure = structure
        self.paired_bases = paired_bases if paired_bases is not None else self.structure_to_paired_bases(structure) if structure is not None else None
        self.dms = dms
        self.opt_dict = {'paired_bases': self.paired_bases, 'dms': self.dms}

    def to_dict(self):
        return self.reference, {'sequence': self.sequence, **{k:v for k,v in self.opt_dict.items() if v is not None}}

    def to_flat_dict(self):
        return {'reference': self.reference, 'sequence': self.sequence, **{k:v for k,v in self.opt_dict.items() if v is not None}}

    def __str__(self):
        return '"'+self.reference+'"' + ':' + str({"sequence": self.sequence, **{k:v for k,v in self.opt_dict.items() if v is not None}}).replace("'",'"').replace('(','[').replace(')',']').replace('None','null')

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.reference}', sequence='{self.sequence}', paired_bases={self.paired_bases}, dms={self.dms})"

    def structure_to_paired_bases(self, structure):
        """Returns a list of tuples (i,j) where i and j are paired bases."""
        paired_bases = []
        stack = []
        for i, char in enumerate(structure):
            if char == '(':
                stack.append(i)
            elif char == ')':
                paired_bases.append((stack.pop(), i))
        return paired_bases

    def embed_sequence(self):
        """Returns a list of integers corresponding to the sequence.

        >>> datapoint = Datapoint(reference='reference', sequence='AACCGG', structure='((..))', dms=[1.0, 2.0, 3.0])
        >>> datapoint.embed_sequence()
        array([1, 1, 2, 2, 3, 3])
        """
        return np.array([seq2int[base] for base in self.sequence])

    def embed_structure(self):
        """Returns a list of integers corresponding to the structure.

        >>> datapoint = Datapoint(reference='reference', sequence='AACCGG', structure='((..))', dms=[1.0, 2.0, 3.0])
        >>> datapoint.embed_structure()
        array([2, 2, 1, 1, 3, 3])
        """
        return np.array([dot2int[base] for base in self.structure])



class DatapointFactory:
    """A datapoint is a tuple (sequence, reference, structure, dms) where:

    - reference is the name of the sequence in the fasta file
    - sequence is a string of A, C, G, U
    - structure is a string of (, ), .
    - dms is a list of floats corresponding to the reactivity of each base in the sequence

    """

    def __call__(sequence, reference, structure=None, dms=None, predict_structure=False, predict_dms=False):
        if structure is None and predict_structure:
            structure = Fasta.predict_structure(sequence)
        return Datapoint(sequence, reference, structure, dms)

    def from_ct(ct_file, predict_dms):
        """Create a datapoint from a ct file. If predict_dms is True, the dms will be predicted using RNAstructure"""
        reference, sequence, paired_bases = Ct.parse(ct_file)
        return Datapoint(sequence, reference, paired_bases=paired_bases, dms = Ct.predict_dms(ct_file) if predict_dms else None)

    def from_fasta(sequence, reference, predict_structure, predict_dms):
        """Create a datapoint from a fasta file. The structure and dms will be None."""
        structure, dms = None, None
        if predict_structure:
            structure = Fasta.predict_structure(sequence)
        if predict_dms:
            dms = Fasta.predict_dms(sequence)
        return Datapoint(sequence, reference, structure, dms)

    def from_dreem_output(reference, sequence, mutation_rate, predict_structure):
        """Create a datapoint from a dreem output file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        return Datapoint(
            sequence=sequence,
            reference=reference,
            structure=Fasta.predict_structure(sequence) if predict_structure else None,
            dms=mutation_rate)


    def from_json_line(string):
        """Create a datapoint from a json line. The json line should have the following format:
        "reference": {"sequence": "sequence", "paired_bases": [[1, 2], [3,4]], "dms": [1.0, 2.0, 3.0]}

        Example:
        >>> DatapointFactory.from_json_line('"reference": {"sequence": "sequence"}')
        Datapoint('reference', sequence='sequence', paired_bases=None, dms=None)
        >>> DatapointFactory.from_json_line('"reference": {"sequence": "sequence", "paired_bases": [[1, 2], [3, 4]], "dms": [1.0, 2.0, 3.0]}')
        Datapoint('reference', sequence='sequence', paired_bases=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0])
        """
        string= string.strip()[:-1] if string.strip()[-1] == ',' else string.strip()
        string = add_braces_if_no_braces(string)
        d = json.loads(string)
        ref = list(d.keys())[0]
        d = d[ref]
        return Datapoint(
            reference=ref,
            sequence=d['sequence'],
            paired_bases=d['paired_bases'] if 'paired_bases' in d else None,
            dms=d['dms'] if 'dms' in d else None
        )



