# RNA_data

A repo to manipulate the data for our RNA structure prediction model.


## Installation

1. Download the code by cloning the repo:

```bash

git clone https://github.com/rouskinlab/RNA_data

```

2. Rename `config_template.yaml` to `my_config.yaml` and change:

- `DATA_FOLDER`, which is the default path for storing the datasets.
- `RNASTRUCTURE_PATH`, you must change this to the path of your RNAstructure executable.
- `RNASTRUCTURE_TEMP_FOLDER`, you can change this to the path of your choice.
- `HUGGINGFACE_TOKEN`, get a token for `huggingface.co/rouskinlab` [here](https://huggingface.co/rouskinlab) and paste it there.

