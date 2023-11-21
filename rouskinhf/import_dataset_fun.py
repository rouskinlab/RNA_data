from os.path import exists
from os import remove
from numpy import ndarray
import numpy as np
from .datafolder import DataFolder
import os

def import_dataset(name:str, force_download:bool=False)->ndarray:

    """Finds the dataset with the given name for the given type of data.

    Parameters
    ----------

    name : str
        Name of the dataset to find.
    force_download : bool
        Whether to force download the dataset from HuggingFace Hub. Defaults to False.

    Returns
    -------

    dict: {str: ndarray}
        Dictionary with the following keys: 'references', 'sequences', and 'structure' or 'DMS' depending on the data type.
    """
    
    # Get the data folder
    try:
        assert not force_download, "Force download from HuggingFace Hub"
        datafolder = DataFolder.from_local(name=name)
        source = 'local'
        print("Using local data for: {}".format(name))
    except AssertionError as e:
        try:
            print(e)
            datafolder = DataFolder.from_huggingface(name=name, overwrite=True)
            source = 'huggingface'
            print("Using data from HuggingFace Hub for {}".format(name))
        except AssertionError as e:
            raise AssertionError("Could not find dataset: {} because of error: {}".format(name, e))

    datafolder.generate_npy()

    out = {
        'references': np.load(datafolder.get_references_npy(), allow_pickle=True),
        'sequences': np.load(datafolder.get_sequences_npy(), allow_pickle=True),
        }

    if exists(datafolder.get_base_pairs_npy()):
        out['base_pairs'] = np.load(datafolder.get_base_pairs_npy(), allow_pickle=True)
    if exists(datafolder.get_dms_npy()):
        out['dms'] = np.load(datafolder.get_dms_npy(), allow_pickle=True)
    if exists(datafolder.get_shape_npy()):
        out['shape'] = np.load(datafolder.get_shape_npy(), allow_pickle=True)

    return out
