from .env import DATA_FOLDER
import pandas as pd
from typing import Dict
import os
from os.path import exists
from huggingface_hub import HfApi
from numpy import ndarray
import numpy as np
from .datafolder import DataFolder

def import_dataset(name:str, data:str, force_download:bool=False)->ndarray:

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

    ndarray
        The dataset with the given name for the given type of data.

    Example
    -------

    >>> import_dataset(name='for_testing', data='structure').shape
    (2,)
    >>> import_dataset(name='for_testing', data='DMS').shape
    (2,)
    >>> import_dataset(name='for_testing', data='structure', force_download=True).shape
    (2,)
    >>> import_dataset(name='for_testing', data='DMS', force_download=True).shape
    (2,)

    """

    assert data in ['structure', 'DMS'], "data must be either 'structure' or 'DMS'"

    # Get the data folder
    try:
        assert not force_download, "Force download from HuggingFace Hub"
        datafolder = DataFolder.from_local(name=name)
    except:
        try:
            datafolder = DataFolder.from_huggingface(name=name)
        except:
            raise ValueError("Dataset not found on HuggingFace Hub")

    # Get the dataset
    if data == 'structure':
        if not exists(datafolder.get_base_pairs_npy()) or force_download:
            datafolder.generate_npy()
        return np.load(datafolder.get_base_pairs_npy(), allow_pickle=True)

    elif data == 'DMS':
        if not exists(datafolder.get_dms_npy()) or force_download:
            datafolder.generate_npy()
        return np.load(datafolder.get_dms_npy(), allow_pickle=True)
