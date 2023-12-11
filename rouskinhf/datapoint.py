from .parsers import *
from .util import (
    standardize_sequence,
    sequence_has_regular_characters,
)
import numpy as np
from .rnastructure import RNAstructure_singleton


class Datapoint:

    """A datapoint is a data structure corresponding to a sequence, a reference, and optionally a structure, dms and shape."""

    def __new__(cls, *args, **kwargs):
        # Check if the conditions are met before creating the object
        sequence = kwargs.get("sequence", args[0] if len(args) > 0 else None)
        reference = kwargs.get("reference", args[1] if len(args) > 1 else None)
        dotbracket = kwargs.get("dotbracket", None)
        if (
            sequence is None
            or reference is None
            or len(sequence) == 0
            or len(reference) == 0
        ):
            return None
        
        if dotbracket is not None:
            if not dotbracket.count("(") == dotbracket.count(")"):
                return None
            
        if not sequence_has_regular_characters(sequence):
            return None

        # If conditions are met, create and return the new instance
        instance = super().__new__(cls)
        return instance

    def __init__(
        self,
        sequence,
        reference,
        dotbracket=None,
        dms=None,
        shape=None,
        structure=None,
    ):
        for attr in [sequence, reference]:
            assert isinstance(
                attr, str
            ), f"Expected {attr} to be a string, got {type(attr)} instead."
            assert len(attr) > 0, f"Expected {attr} to be non-empty."

        # standardize the sequence
        sequence = standardize_sequence(sequence)

        if not sequence_has_regular_characters(sequence):
            raise Exception(
                f"Sequence {sequence} contains characters other than ACGTUacgtu."
            )

        self.reference = reference
        self.sequence = sequence
        self.structure = structure
        
        if structure is None and dotbracket is not None:
            self.structure = self.dotbracket_to_structure(dotbracket)
        if dms is not None:
            self.dms = self._format_signal(dms)
        if shape is not None:
            self.shape = self._format_signal(shape)

    def get_opt_dict(self):
        return {
            attr: eval(f"self.{attr}")
            for attr in ["structure", "dms", "shape"]
            if hasattr(self, attr)
        }

    def _assert_structure(self):
        """Ensures that the base pairs are 0=<bp<len(sequence) and that that there's maximum one pair per base.
        >>> datapoint = Datapoint(reference='reference', sequence='AACCGG', dotbracket='((..))')
        >>> datapoint._assert_structure()
        True
        >>> datapoint = Datapoint(reference='reference', sequence='AACCGG', dotbracket='((..)(')
        >>> datapoint = Datapoint(reference='reference', sequence='AACCGG', structure=[[1, 2], [3, 4]])
        >>> datapoint._assert_structure()
        True
        >>> datapoint = Datapoint(reference='reference', sequence='AACCGG', structure=[[1, 2], [3, 3]])
        >>>> datapoint._assert_structure()
        False
        >>> datapoint = Datapoint(reference='reference', sequence='AACCGG', structure=[[1, 2], [3, 6]])
        >>> datapoint._assert_structure()
        False   
        >>> datapoint = Datapoint(reference='reference', sequence='AACCGG', structure=[[1, 2], [3, -1]])
        >>> datapoint._assert_structure()
        False
        """
        sequence = self.sequence
        structure = self.structure

        if structure is None:
            return True
        return (
            all(
                # check boundaries
                [
                    0 <= pair[0] < len(sequence)
                    and 0 <= pair[1] < len(sequence)
                    and pair[0] != pair[1]
                    for pair in structure
                ]
            )
            # check that there's maximum one pair per base
            and len(
                set(
                    [pair[0] for pair in structure]
                    + [pair[1] for pair in structure]
                )
            )
            == len(structure) * 2
        )

    def convert_arrays_to_tuple(self):
        if hasattr(self, "dms") and self.dms is not None:
            self.dms = tuple(self.dms)
        if hasattr(self, "shape") and self.shape is not None:
            self.shape = tuple(self.shape)
        if hasattr(self, "structure") and self.structure is not None:
            self.structure = tuple([tuple(pair) for pair in self.structure])

    def convert_arrays_to_list(self):
        if hasattr(self, "dms") and self.dms is not None:
            self.dms = list(self.dms)
        if hasattr(self, "shape") and self.shape is not None:
            self.shape = list(self.shape)
        if hasattr(self, "structure") and self.structure is not None:
            self.structure = [list(pair) for pair in self.structure]

    def _format_signal(self, signal):
        return [round(d, 4) for d in signal]

    def to_dict(self):
        return self.reference, {
            "sequence": self.sequence,
            **{k: v for k, v in self.get_opt_dict().items() if v is not None},
        }

    def to_flat_dict(self):
        return {
            "reference": self.reference,
            "sequence": self.sequence,
            **{k: v for k, v in self.get_opt_dict().items() if v is not None},
        }

    @classmethod
    def from_flat_dict(cls, d):
        return cls(
            sequence=d["sequence"],
            reference=d["reference"],
            **{k: v for k, v in d.items() if k not in ["sequence", "reference"]},
        )

    def __str__(self):
        return f'"{self.reference}":' + str(
            {
                "sequence": self.sequence,
                **{k: v for k, v in self.get_opt_dict().items() if v is not None},
            }
        ).replace("'", '"').replace("(", "[").replace(")", "]").replace("None", "null")

    def __repr__(self) -> str:
        out = (
            f"{self.__class__.__name__}('{self.reference}', sequence='{self.sequence}'"
        )
        if hasattr(self, "structure"):
            out += f", structure={self.structure}"
        if hasattr(self, "dms"):
            out += f", dms={self.dms}"
        return out + ")"

    def dotbracket_to_structure(self, dotbracket):
        """Returns a list of tuples (i,j) where i and j are paired bases."""
        structure = set()
        stack = []
        for i, char in enumerate(dotbracket):
            if char == "(":
                stack.append(i)
            elif char == ")":
                structure.add((stack.pop(), i))
        return structure


class DatapointFactory:
    """Factory class to create datapoints from different formats."""

    def from_bpseq(bpseq_file):
        """Create a datapoint from a bpseq file. If predict_dms is True, the dms will be predicted using RNAstructure"""
        reference, sequence, structure = BPseq.parse(bpseq_file)
        sequence = standardize_sequence(sequence)

        if sequence_has_regular_characters(sequence):
            return Datapoint(
                sequence,
                reference,
                structure=structure,
            )

    def from_ct(ct_file):
        """Create a datapoint from a ct file. If predict_dms is True, the dms will be predicted using RNAstructure"""
        reference, sequence, structure = Ct.parse(ct_file)
        sequence = standardize_sequence(sequence)

        if sequence_has_regular_characters(sequence):
            return Datapoint(
                sequence,
                reference,
                structure=structure,
            )

    def from_fasta(sequence, reference, predict_structure):
        """Create a datapoint from a fasta file. The structure and dms will be None."""
        sequence = standardize_sequence(sequence)

        if sequence_has_regular_characters(sequence):
            dotbracket, dms = None, None
            if predict_structure:
                dotbracket = Fasta.predict_structure(sequence)
            return Datapoint(sequence, reference, dotbracket=dotbracket, dms=dms)

    def from_dreem_output(reference, sequence, mutation_rate, predict_structure):
        """Create a datapoint from a dreem output file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        sequence = standardize_sequence(sequence)
        mutation_rate = np.array([float(m) for m in mutation_rate], dtype=np.float32)
        if sequence_has_regular_characters(sequence):
            return Datapoint(
                sequence=sequence,
                reference=reference,
                dotbracket=RNAstructure_singleton.predictStructure(
                    sequence, dms=mutation_rate
                )
                if predict_structure
                else None,
                dms=mutation_rate,
            )

    def from_json_line(ref, d, predict_structure=False):
        """Create a datapoint from a json line. The json line should have the following format:
        "reference": {"sequence": "sequence", "structure": [[1, 2], [3,4]], "dms": [1.0, 2.0, 3.0]}

        Example:
        >>> DatapointFactory.from_json_line("reference", {"sequence": "ACAAGU"})
        Datapoint('reference', sequence='ACAAGU')
        >>> DatapointFactory.from_json_line("reference", {"sequence": "ACAAGU", "structure": [[1, 2], [3, 4]], "dms": [1.0, 2.0, 3.0]})
        Datapoint('reference', sequence='ACAAGU', structure=((1, 2), (3, 4)), dms=(1.0, 2.0, 3.0))
        >>> DatapointFactory.from_json_line("reference", {"sequence": "something else than ACGTUacgtu", "structure": [[1, 2], [3, 4]], "dms": [1.0, 2.0, 3.0]})
        """

        # create the datapoint
        sequence = d["sequence"]
        sequence = standardize_sequence(sequence)

        # hack
        if predict_structure and (not "structure" in d) and (not "dotbracket" in d):
            d["dotbracket"] = RNAstructure_singleton.predictStructure(
                sequence, dms=d["dms"] if "dms" in d else None
            )

        if sequence_has_regular_characters(sequence):
            return Datapoint(
                reference=ref,
                sequence=sequence,
                structure=d["structure"] if "structure" in d else None,
                dotbracket=d["dotbracket"] if "dotbracket" in d else None,
                dms=d["dms"] if "dms" in d else None,
                shape = d["shape"] if "shape" in d else None,
            )
