[![Python CI](https://github.com/rouskinlab/rouskinhf/actions/workflows/CI.yml/badge.svg)](https://github.com/rouskinlab/rouskinhf/actions/workflows/CI.yml)
[![Publish distributions 📦 to PyPI](https://github.com/rouskinlab/rouskinhf/actions/workflows/release.yml/badge.svg)](https://github.com/rouskinlab/rouskinhf/actions/workflows/release.yml)
![PyPI](https://img.shields.io/pypi/v/rouskinhf)
![GitHub tag (with filter)](https://img.shields.io/github/v/tag/rouskinlab/rouskinhf)

# Download your RNA data from HuggingFace with rouskinhf!

A repo to manipulate the data for our RNA structure prediction model. This repo allows you to:
- pull datasets from the Rouskinlab's HuggingFace
- create datasets from local files and push them to HuggingFace, from the formats:
    - `.fasta`
    - `.ct`
    - `.json` (DREEM output format)
    - `.json` (Rouskinlab's huggingface format)

## Important notes

- Sequences with bases different than `A`, `C`, `G`, `T`, `U`, `N`, `a`, `c`, `g`, `t`, `u`, `n` are not supported. The data will be filtered out.

## Dependencies
- [RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html) (also available on [Rouskinlab GitHub](https://github.com/rouskinlab/RNAstructure)).

## Push a new release to Pypi

1. Edit version to `vx.y.z` in `pyproject.toml`. Then run in a terminal `git add . && git commit -m 'vx.y.z' && git push`.
2. Create and push a git tag `vx.y.z` by running in a terminal `git tag 'vx.y.z' && git push --tag`.
3. Create a release for the tag `vx.y.z` on Github Release.
4. Make sure that the Github Action `Publish distributions 📦 to PyPI` passed on Github Actions.

## Installation

### Get a HuggingFace token

Go to [HuggingFace](https://huggingface.co/) and create an account. Then go to your profile and copy your token ([huggingface.co/settings/tokens](https://huggingface.co/settings/tokens)).

### Create an environment file

Open a terminal and type:

```bash
nano env
```

Copy paste the following content, and change the values to your own:

```bash
export HUGGINGFACE_TOKEN="your token here"  # you must change this to your HuggingFace token
export DATA_FOLDER="data/datafolders" # where the datafolder are stored by default, change it if you want to store it somewhere else
export DATA_FOLDER_TESTING="data/input_files_for_testing" # Don't touch this
export RNASTRUCTURE_PATH="/Users/ymdt/src/RNAstructure/exe" # Change this to the path of your RNAstructure executable
export RNASTRUCTURE_TEMP_FOLDER="temp" # You can change this to the path of your RNAstructure temp folder
```

Then save the file and exit nano.

### Source the environment

```bash
source env
```

### Install the package with pip

```bash
pip install rouskinhf
```


## Tutorials

### Authentify your machine to HuggingFace

See the [tutorial](https://github.com/rouskinlab/rouskinhf/blob/main/tutorials/huggingface.ipynb).

### Download a datafolder from HuggingFace

See the [tutorial](https://github.com/rouskinlab/rouskinhf/blob/main/tutorials/use_for_models.ipynb).

### Create a datafolder from local files and push it to HuggingFace

See the [tutorial](https://github.com/rouskinlab/rouskinhf/blob/main/tutorials/create_push_pull.ipynb).

## About

### Sourcing the environment and keeping your environment variable secret

The variables defined in the `env` file are required by `rouskinhf`. Make that before you use `rouskinhf`, you run in a terminal:

```bash
source env
```
 or, in a Jupyter notebook:

```python
!pip install python-dotenv
%load_ext dotenv
%dotenv env
```

 The point of using environment variables is to ensure the privacy of your huggingface token. Make sure to add your `env` file to your `.gitignore`, so your HuggingFace token doesn't get pushed to any public repository.

### Import data with ``import_dataset``

This repo provides a function ``import_dataset``, which allows your to pull a dataset from HuggingFace and store it locally. If the data is already stored locally, it will be loaded from the local folder. The type of data available is the DMS signal and the structure, under the shape of paired bases tuples. The function has the following signature:

```python
def import_dataset(name:str, data:str, force_download:bool=False)->np.ndarray:

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

    >>> import_dataset(name='for_testing', data='structure').keys()
    dict_keys(['references', 'sequences', 'structure'])
    >>> import_dataset(name='for_testing', data='DMS').keys()
    dict_keys(['references', 'sequences', 'DMS'])
    >>> import_dataset(name='for_testing', data='structure', force_download=True).keys()
    dict_keys(['references', 'sequences', 'structure'])
    >>> import_dataset(name='for_testing', data='DMS', force_download=True).keys()
    dict_keys(['references', 'sequences', 'DMS'])
```

### FYI, the datafolder object

The datafolder object is a wrapper around your local folder and HuggingFace API, to keep a consistent datastructure across your datasets. It contains multiple methods to create datasets from various input formats, store the data and metadata in a systematic way, and push / pull from HuggingFace.

On HuggingFace, the datafolder stores the data under the following structure:

```bash
HUGGINGFACE DATAFOLDER
- [datafolder name]
    - source
        - whichever file(s) you used to create the dataset (fasta, set of CTs, etc.).
    - data.json # the data under a human readable format.
    - info.json # the metadata of the dataset. This file indicates how we got the DMS signal and the structures (directly from the source or from a prediction).
    - README.md # the metadata of the dataset in a human readable format.
```

Locally, we have the same structure with the addition of .npy files which contain the data in a machine readable format. Each .npy file contains a numpy array of the data, and the name of the file is the name of the corresponding key in the data.json file. The source file won’t be downloaded by default. Hence, the local structure is:

```bash
LOCAL DATAFOLDER
- [datafolder name]
    ...
    - README.md # the metadata of the dataset in a human readable format
    - references.npy
    - sequences.npy
    - base_pairs.npy
    - dms.npy
```
