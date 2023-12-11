import os
import numpy as np
from .datapoint import Datapoint, DatapointFactory
from typing import List
from .parsers import Fasta, DreemOutput

import pandas as pd
from tqdm import tqdm as tqdm_parser
import json


class ListofDatapoints:
    """Class to store a list of datapoints."""
    def __init__(self, datapoints=[], verbose=True):
        self.datapoints = datapoints

    def __call__(self) -> List[Datapoint]:
        return self.datapoints
    
    def __len__(self) -> int:
        return len(self.datapoints)

    @classmethod
    def from_fasta(
        cls, fasta_file, predict_structure, tqdm=True, verbose=True
    ) -> "ListofDatapoints":
        """Create a list of datapoint from a fasta file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        sequences, references = Fasta.parse(fasta_file)
        return cls(
            [
                DatapointFactory.from_fasta(sequence, reference, predict_structure)
                for sequence, reference in tqdm_parser(
                    zip(sequences, references),
                    total=len(sequences),
                    desc="Parsing fasta file",
                    disable=not tqdm,
                )
            ],
            verbose=verbose,
        )

    @classmethod
    def from_bpseq(cls, bpseq_files, tqdm=True, verbose=True):
        """Create a list of datapoint from a bpseq file. The dms will be predicted if predict_dms is True."""
        return cls(
            [
                DatapointFactory.from_bpseq(ct_file)
                for ct_file in tqdm_parser(
                    bpseq_files,
                    total=len(bpseq_files),
                    desc="Parsing bpseq files",
                    disable=not tqdm,
                )
            ],
            verbose=verbose,
        )

    @classmethod
    def from_ct(cls, ct_files, tqdm=True, verbose=True):
        """Create a list of datapoint from a list of ct files. The dms will be predicted if predict_dms is True."""
        return cls(
            [
                DatapointFactory.from_ct(ct_file)
                for ct_file in tqdm_parser(
                    ct_files,
                    total=len(ct_files),
                    desc="Parsing ct files",
                    disable=not tqdm,
                )
            ],
            verbose=verbose,
        )

    @classmethod
    def from_dreem_output(
        cls, dreem_output_file, predict_structure, tqdm=True, verbose=True
    ):
        """Create a list of datapoint from a dreem output file. The structure and dms will be predicted if predict_structure and predict_dms are True."""
        n_lines = len(list(DreemOutput.parse(dreem_output_file)))
        return cls(
            [
                DatapointFactory.from_dreem_output(
                    reference, sequence, mutation_rate, predict_structure
                )
                for reference, sequence, mutation_rate in tqdm_parser(
                    DreemOutput.parse(dreem_output_file),
                    total=n_lines,
                    desc="Parsing dreem output file",
                    disable=not tqdm,
                )
            ],
            verbose=verbose,
        )

    @classmethod
    def from_json(cls, json_file, predict_structure=False, tqdm=True, verbose=True):
        """Create a list of datapoint from a json file."""
        data = json.load(open(json_file))
        return cls(
            [
                DatapointFactory.from_json_line(reference, line, predict_structure)
                for reference, line in tqdm_parser(
                    data.items(),
                    total=len(data),
                    desc="Parsing json file",
                    disable=not tqdm,
                )
            ],
            verbose=verbose,
        )

    def to_dict(self) -> dict:
        """Converts the list of datapoints into a dictionary."""
        return {
            datapoint.reference: datapoint.to_dict()[1] for datapoint in self.datapoints
        }

    def to_json(self, path) -> None:
        """Write a list of datapoints to a json file."""
        if os.path.dirname(path) != "":
            os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write("{\n")
            for idx, datapoint in enumerate(
                self.datapoints
            ):  # write the datapoints one by one to avoid memory issues
                f.write(
                    str(datapoint) + ",\n"
                    if idx != len(self.datapoints) - 1
                    else str(datapoint) + "\n"
                )
            f.write("}\n")

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
        if datapoints is None:
            datapoints = self.datapoints
        for dp in datapoints:
            dp.convert_arrays_to_tuple()
        return pd.DataFrame([datapoint.to_flat_dict() for datapoint in datapoints])

    def from_pandas(self, df: pd.DataFrame) -> None:
        """Converts a pandas dataframe into a list of datapoints.

        Example:
            >>> df = pd.DataFrame([{'reference': 'reference', 'sequence': 'AACCGG', 'structure': [[1, 2], [3, 4]], 'dms': [1.0, 2.0, 3.0]}])
            >>> ListofDatapoints(verbose=False).from_pandas(df)
            [Datapoint('reference', sequence='AACCGG', structure=[[1, 2], [3, 4]], dms=[1.0, 2.0, 3.0])]
        """
        datapoints = [
            Datapoint.from_flat_dict(datapoint_dict)
            for datapoint_dict in df.to_dict("records")
        ]
        for dp in datapoints:
            dp.convert_arrays_to_list()
        return datapoints

    def drop_none_dp(self):
        n_input_datapoints = len(self.datapoints)
        self.datapoints = list(filter(lambda x: x is not None, self.datapoints))
        n_unvalid_datapoints = n_input_datapoints - len(self.datapoints)
        return self.datapoints, n_unvalid_datapoints

    
