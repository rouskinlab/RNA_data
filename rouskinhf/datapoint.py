
from typing import Any, List
from .parsers import *
from .util import add_braces_if_no_braces, dot2int, int2dot, seq2int, standardize_sequence, sequence_has_regular_characters
import numpy as np
from .rnastructure import RNAstructure_singleton

class Datapoint:

    """A datapoint is a data structure where:
    - reference is the name of the sequence in the fasta file
    - sequence is a string of A, C, G, U
    - structure is a string of (, ), .
    - dms is a list of floats corresponding to the reactivity of each base in the sequence
    - paired_bases is a list of tuples (i,j) where i and j are paired bases.


    Example:
    >>> print(Datapoint(reference='reference', sequence='AACCGG', structure='((..))', dms=[1.0, 2.0, 3.0]))
    "reference":{"sequence": "AACCGG", "paired_bases": [[0, 5], [1, 4]], "dms": [1.0, 2.0, 3.0]}
    >>> print(Datapoint(reference='reference', sequence='not a valid sequence', structure='((..))'))
    None
    >>> print(Datapoint(reference='reference', sequence='NNNNNNN', structure='((..))'))
    None
    >>> print(Datapoint(reference='reference', sequence='AACCGG', paired_bases=[(0, 5), (1, 4)]))
    "reference":{"sequence": "AACCGG", "paired_bases": [[0, 5], [1, 4]]}
    >>> print(Datapoint(reference='reference', sequence='AACCGG', paired_bases=[(1, 5), (1, 4)]))
    None
    """

    def __new__(cls, *args, **kwargs):
        # Check if the conditions are met before creating the object
        sequence = kwargs.get('sequence', args[0] if len(args) > 0 else None)
        reference = kwargs.get('reference', args[1] if len(args) > 1 else None)
        if sequence is None or reference is None:
            return None
        sequence = standardize_sequence(sequence)

        if not sequence_has_regular_characters(sequence) or len(sequence) == 0 or len(reference) == 0:
            return None
        
        if 'paired_bases' in kwargs:
            if not cls._assert_paired_bases(None, kwargs['paired_bases'], sequence):
                return None

        # If conditions are met, create and return the new instance
        instance = super().__new__(cls)
        return instance

    def __init__(self, sequence, reference, structure=None, dms=None, paired_bases=None):
        for attr in [sequence, reference]:
            assert isinstance(attr, str), f"Expected {attr} to be a string, got {type(attr)} instead."
            assert len(attr) > 0, f"Expected {attr} to be non-empty."

        # standardize the sequence
        sequence = standardize_sequence(sequence)

        if not sequence_has_regular_characters(sequence):
            raise Exception(f"Sequence {sequence} contains characters other than ACGTUacgtu.")

        self.reference = reference
        self.sequence = sequence
        self.structure = structure
        if paired_bases is not None or structure is not None:
            self.paired_bases = self._format_paired_bases(paired_bases) if paired_bases is not None else self._format_paired_bases(self.structure_to_paired_bases(structure))
        if dms is not None:
            self.dms = self._format_dms(dms) 
        self.opt_dict = {attr: eval(f'self.{attr}') for attr in ['paired_bases', 'dms'] if hasattr(self, attr)}

    def _format_paired_bases(self, paired_bases):
        """Returns a set of tuples.
        >>> datapoint = Datapoint(reference='reference', sequence='AACCGG', structure='((..))', dms=[1.0, 2.0, 3.0])
        >>> datapoint._format_paired_bases([[1, 4], [0, 5]])
        ((0, 5), (1, 4))
        """
        return tuple(sorted([tuple(sorted(pair)) for pair in paired_bases])) if paired_bases is not None else None
    
    def _assert_paired_bases(self, paired_bases, sequence):
        """Ensures that the base pairs are 0=<bp<len(sequence) and that that there's maximum one pair per base. 
        >>> datapoint = Datapoint(reference='reference', sequence='AACCGG', structure='((..))', dms=[1.0, 2.0, 3.0])
        >>> datapoint._assert_paired_bases([(0, 5), (1, 4)], 'AACCGG')
        True
        >>> datapoint._assert_paired_bases([(0, 5), (1, 4), (1, 2)], 'AACCGG')
        False
        >>> datapoint._assert_paired_bases([(0, 5), (4, 1), (1, 2)], 'AACCGG')
        False
        >>> datapoint._assert_paired_bases([(0, 8), (3, 4), (1, 2)], 'AACCGG')
        False
        >>> datapoint._assert_paired_bases(None, 'AACCGG')
        True
        """
        if paired_bases is not None:
            return all([0 <= pair[0] < len(sequence) and 0 <= pair[1] < len(sequence) and pair[0] != pair[1] for pair in paired_bases]) \
                and len(set([pair[0] for pair in paired_bases] + [pair[1] for pair in paired_bases])) == len(paired_bases) * 2
        else:
            return True

    def _format_dms(self, dms):
        """returns a tuple
        >>> datapoint = Datapoint(reference='reference', sequence='AACCGG', structure='((..))', dms=[1.0, 2.0, 3.0])
        >>> datapoint._format_dms([1.0, 2.0, 3.0])
        (1.0, 2.0, 3.0)
        """
        if type(dms) == np.ndarray:
            dms = dms.tolist()
        return tuple(dms) 
        

    def to_dict(self):
        return self.reference, {'sequence': self.sequence, **{k:v for k,v in self.opt_dict.items() if v is not None}}

    def to_flat_dict(self):
        return {'reference': self.reference, 'sequence': self.sequence, **{k:v for k,v in self.opt_dict.items() if v is not None}}

    @classmethod
    def from_flat_dict(cls, d):
        return cls(sequence=d['sequence'], reference=d['reference'], **{k:v for k,v in d.items() if k not in ['sequence', 'reference']})

    def __str__(self):
        return '"'+self.reference+'"' + ':' + str({"sequence": self.sequence, **{k:v for k,v in self.opt_dict.items() if v is not None}}).replace("'",'"').replace('(','[').replace(')',']').replace('None','null')

    def __repr__(self) -> str:
        out = f"{self.__class__.__name__}('{self.reference}', sequence='{self.sequence}'"
        if hasattr(self, 'paired_bases'):
            out += f", paired_bases={self.paired_bases}"
        if hasattr(self, 'dms'):
            out += f", dms={self.dms}"
        return out + ')'

    def structure_to_paired_bases(self, structure):
        """Returns a list of tuples (i,j) where i and j are paired bases."""
        paired_bases = set()
        stack = []
        for i, char in enumerate(structure):
            if char == '(':
                stack.append(i)
            elif char == ')':
                paired_bases.add((stack.pop(), i))
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

        >>> from numpy import array
        >>> datapoint = Datapoint(reference='reference', sequence='AACCGG', structure='((..))', dms=[1.0, 2.0, 3.0])
        >>> assert not (datapoint.embed_structure() - array([2, 2, 1, 1, 3, 3])).any(), 'The sequence is not embedded correctly.'
        """
        return np.array([dot2int[base] for base in self.structure])



class DatapointFactory:
    """A datapoint is a tuple (sequence, reference, structure, dms) where:

    - reference is the name of the sequence in the fasta file
    - sequence is a string of A, C, G, U
    - structure is a string of (, ), .
    - dms is a list of floats corresponding to the reactivity of each base in the sequence

    """
    def from_ct(ct_file, predict_dms):
        """Create a datapoint from a ct file. If predict_dms is True, the dms will be predicted using RNAstructure"""
        reference, sequence, paired_bases = Ct.parse(ct_file)
        sequence = standardize_sequence(sequence)

        if sequence_has_regular_characters(sequence):
            return Datapoint(sequence, reference, paired_bases=paired_bases, dms = Ct.predict_dms(ct_file) if predict_dms else None)

    def from_fasta(sequence, reference, predict_structure, predict_dms):
        """Create a datapoint from a fasta file. The structure and dms will be None."""
        sequence = standardize_sequence(sequence)

        if sequence_has_regular_characters(sequence):
            structure, dms = None, None
            if predict_structure:
                structure = Fasta.predict_structure(sequence)
            if predict_dms:
                dms = Fasta.predict_dms(sequence)
            return Datapoint(sequence, reference, structure, dms)

    def from_dreem_output(reference, sequence, mutation_rate, predict_structure):
        """Create a datapoint from a dreem output file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        sequence = standardize_sequence(sequence)

        if sequence_has_regular_characters(sequence):
            return Datapoint(
                sequence=sequence,
                reference=reference,
                structure=RNAstructure_singleton.predictStructure(sequence, dms=mutation_rate) if predict_structure else None,
                dms=mutation_rate)


    def from_json_line(string):
        """Create a datapoint from a json line. The json line should have the following format:
        "reference": {"sequence": "sequence", "paired_bases": [[1, 2], [3,4]], "dms": [1.0, 2.0, 3.0]}

        Example:
        >>> DatapointFactory.from_json_line('"reference": {"sequence": "ACAAGU"}')
        Datapoint('reference', sequence='ACAAGU')
        >>> DatapointFactory.from_json_line('"reference": {"sequence": "ACAAGU", "paired_bases": [[1, 2], [3, 4]], "dms": [1.0, 2.0, 3.0]}')
        Datapoint('reference', sequence='ACAAGU', paired_bases=((1, 2), (3, 4)), dms=(1.0, 2.0, 3.0))
        >>> DatapointFactory.from_json_line('"reference": {"sequence": "something else than ACGTUacgtu", "paired_bases": [[1, 2], [3, 4]], "dms": [1.0, 2.0, 3.0]}')
        """

        # process the string into a dictionary
        string= string.strip()[:-1] if string.strip()[-1] == ',' else string.strip()
        string = add_braces_if_no_braces(string)
        d = json.loads(string)
        ref = list(d.keys())[0]
        d = d[ref]

        # create the datapoint
        sequence = d['sequence']
        sequence = standardize_sequence(sequence)

        if sequence_has_regular_characters(sequence):
            return Datapoint(
                reference=ref,
                sequence=sequence,
                paired_bases=d['paired_bases'] if 'paired_bases' in d else None,
                dms=d['dms'] if 'dms' in d else None
            )



