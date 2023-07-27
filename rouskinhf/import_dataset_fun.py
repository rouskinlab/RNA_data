from os.path import exists
from numpy import ndarray
import numpy as np
from .datafolder import DataFolder

def import_dataset(name:str, data:str, force_download:bool=False, force_generate_npy=False)->ndarray:

    """Finds the dataset with the given name for the given type of data.

    Parameters
    ----------

    name : str
        Name of the dataset to find.
    data : str
        Name of the type of data to find the dataset for (structure or DMS).
    force_download : bool
        Whether to force download the dataset from HuggingFace Hub. Defaults to False.

    Returns
    -------

    dict: {str: ndarray}
        Dictionary with the following keys: 'references', 'sequences', and 'structure' or 'DMS' depending on the data type.

    Example
    -------

    >>> import_dataset(name='for_testing', data='structure', force_download=True).keys()
    Force download from HuggingFace Hub
    Over a total of 2 datapoints, there are:
        - 2 valid datapoints
        - 0 invalid datapoints (ex: sequence with non-regular characters)
        - 0 datapoints with the same reference
        - 0 duplicate sequences with the same structure / dms
        - 0 duplicate sequences with different structure / dms
    Using data from HuggingFace Hub
    dict_keys(['references', 'sequences', 'structure'])
    """
    if data == 'dms': data = 'DMS'
    assert data in ['structure', 'DMS'], "data must be either 'structure' or 'DMS'"

    # Get the data folder
    try:
        assert not force_download, "Force download from HuggingFace Hub"
        datafolder = DataFolder.from_local(name=name, generate_npy=force_generate_npy)
        source = 'local'
        print("Using local data")
    except AssertionError as e:
        try:
            print(e)
            datafolder = DataFolder.from_huggingface(name=name)
            source = 'huggingface'
            print("Using data from HuggingFace Hub")
        except:
            raise ValueError("Dataset not found on HuggingFace Hub")

    if not exists(datafolder.get_references_npy()) or not exists(datafolder.get_sequences_npy()) or source == 'huggingface' or force_generate_npy:
        datafolder.generate_npy()

    out = {
        'references': np.load(datafolder.get_references_npy(), allow_pickle=True),
        'sequences': np.load(datafolder.get_sequences_npy(), allow_pickle=True)
        }

    # Get the dataset
    if data == 'structure':
        if not exists(datafolder.get_base_pairs_npy()) or force_download:
            datafolder.generate_npy()
        out['structure'] = np.load(datafolder.get_base_pairs_npy(), allow_pickle=True)

    elif data == 'DMS':
        if not exists(datafolder.get_dms_npy()) or force_download:
            datafolder.generate_npy()
        out['DMS'] = np.load(datafolder.get_dms_npy(), allow_pickle=True)

    return out
