![PyPI](https://img.shields.io/pypi/v/rouskinhf)
![GitHub tag (with filter)](https://img.shields.io/github/v/tag/rouskinlab/rouskinhf)

# Download your RNA data from HuggingFace with rouskinhf!

A wrapper around Huggingface the load data for eFold. You can:
- pull datasets from the Rouskinlab's HuggingFace
- create datasets from local files 

# Installation

### To download data

```bash
pip install rouskinhf
```

### To push data to huggingface (optional) 

- get a token access from the rouskilab huggingface's page
- add this token to your environment

```bash
export HUGGINGFACE_TOKEN="hf_yourtokenhere"
```

### To predict structures from rouskinhf (optional)
You'll need to install D. Mathew's [RNAstructure Fold](https://rna.urmc.rochester.edu/RNAstructure.html) (also available on [Rouskinlab GitHub](https://github.com/rouskinlab/RNAstructure)).

Check your RNAstructure Fold installation in a terminal:

```bash
Fold --version
```

# How to use

### Download a dataset

```python
import rouskinhf

rouskinhf.get_dataset(
    name='bpRNA-1m', # the name of a dataset from huggingface/rouskinlab
    force_download = False # use a local copy of the data if it exists
)
```

### Convert whatever format to rouskinhf format

```python
import rouskinhf

rouskinhf.convert(
    format = 'ct', # can be ct, seismic, bpseq, fasta or json (rouskinhf output data structure)
    file_or_folder = 'path/to/my/ct/folder',
    predict_structure = False, # Add structure from RNAstructure
    filter = True, # removes duplicates, non-regular characters and low AUROC
    min_AUROC=0.8,
)
```
> Note: Sequences with bases different than `A`, `C`, `G`, `T`, `U`, `N`, `a`, `c`, `g`, `t`, `u`, `n` are not supported. The data will be filtered out.


### Rouskinhf structure format
```json
# rouskinhf_output_file.json
{
    "reference_name": {
        "sequence": "CACGCUAUG",
        "structure": [(0,8), (1,7)], # base pair representation
        # whatever other info you need
    }
}
```


