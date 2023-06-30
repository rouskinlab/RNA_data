
from typing import Any, List
from .file_handlers import *

class Datapoint:

    def __init__(self, sequence, reference, structure=None, dms=None, paired_bases=None):
        self.reference = reference
        self.sequence = sequence
        self.paired_bases = paired_bases if paired_bases is not None else self.structure_to_paired_bases(structure) if structure is not None else None
        self.dms = dms
        self.opt_dict = {'paired_bases': self.paired_bases, 'dms': self.dms}

    def to_dict(self):
        return self.reference, {'sequence': self.sequence, **{k:v for k,v in self.opt_dict.items() if v is not None}}

    def __str__(self):
        # TODO this could be improved
        return '"'+self.reference+'"' + ':' + str({"sequence": self.sequence, **{k:v for k,v in self.opt_dict.items() if v is not None}}).replace("'",'"').replace('(','[').replace(')',']').replace('None','null')

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


class ListofDatapoints:

    def from_fasta(fasta_file, predict_structure, predict_dms):
        """Create a list of datapoint from a fasta file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        sequences, references = Fasta.parse(fasta_file)
        return [
            DatapointFactory.from_fasta(sequence, reference, predict_structure, predict_dms)
            for sequence, reference in zip(sequences, references)
        ]

    def from_ct(ct_files, predict_dms):
        """Create a list of datapoint from a ct file. The dms will be predicted if predict_dms is True."""
        return [DatapointFactory.from_ct(ct_file, predict_dms) for ct_file in ct_files]

    def from_dreem_output(dreem_output_file, predict_structure):
        """Create a list of datapoint from a dreem output file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        return [
            DatapointFactory.from_dreem_output(reference, sequence, mutation_rate, predict_structure)
            for reference, sequence, mutation_rate in DreemOutput.parse(dreem_output_file)
        ]


def write_list_of_datapoints_to_json(path, datapoints:List[Datapoint]) -> None:
    """Write a list of datapoints to a json file."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write('{\n')
        for idx, datapoint in enumerate(datapoints): # write the datapoints one by one to avoid memory issues
            f.write(str(datapoint)+ ',\n' if idx != len(datapoints)-1 else str(datapoint)+'\n')
        f.write('}\n')

