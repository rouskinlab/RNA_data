from .env import DATA_FOLDER
import pandas as pd
from typing import Dict
import os
from huggingface_hub import HfApi


class Finder:
    """This class is used to find the datafolders locally and on Hugging Face."""

    def to_pandas(root=DATA_FOLDER)->pd.DataFrame:
        """Returns the list of datafolders."""
        local_datafolders = [d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))]

        api = HfApi()


        #raise NotImplementedError

    def filter_by(where='all', structure=None, dms=None, root=DATA_FOLDER)->Dict[str, list]:
        """Returns a new Finder object with datapoints filtered by the given parameters.

        Args:
            where (str): where the datafolder is stored. Can be 'local', 'huggingface', or 'all'. Defaults to 'all'.

            structure (str): where the RNA structures come from. Can be either 'source', 'RNAstructure', 'all', or 'None'.
                "source" means that the RNA structures come from the source data.
                "RNAstructure" means that the RNA structures come from RNAstructure. If the DMS is 'source', RNAstructure used the DMS signal for the prediction.
                "all" means that the RNA structures come from both the source data and RNAstructure.
                None means that we don't filter by structure source. This is the default.

            dms (str): where the DMS signal comes from. Can be either 'source', 'RNAstructure', 'all', or 'None'.
                "source" means that the DMS signal comes from the source data.
                "RNAstructure" means that the DMS signal comes from RNAstructure.
                "all" means that the DMS signal comes from both the source data and RNAstructure.
                None means that we don't filter by DMS source. This is the default.

            root (str): the root directory to search for datafolders. Defaults to the DATA_FOLDER environment variable.

        Returns:
            A dictionary of datafolders names:
                {
                    'local': ['datafolder1', 'datafolder2', ...],
                    'huggingface': ['datafolder1', 'datafolder2', ...]
                }
        """

        df = Finder.to_pandas(root=root)

        assert where in ['local', 'huggingface', 'all'], "where must be either 'local', 'huggingface', or 'all'."
        assert structure in ['source', 'RNAstructure', 'all', None], "structure must be either 'source', 'RNAstructure', 'all', or None."
        assert dms in ['source', 'RNAstructure', 'all', None], "dms must be either 'source', 'RNAstructure', 'all', or None."

        #raise NotImplementedError
        return {'local': [], 'huggingface': []}