from huggingface_hub import snapshot_download
import os
from os.path import dirname, exists
import json
from huggingface_hub import HfApi
import datetime
import numpy as np
import json

from .path import Path
from .env import Env


def get_dataset(name: str, path='data', force_download=False, tqdm=True):
    """Get a dataset from HuggingFace or from the local cache.

    Args:
        name (str): Name of the dataset.
        force_download (bool, optional): Whether to force the download or not. Defaults to False.
        tqdm (bool, optional): Whether to display a progress bar or not. Defaults to True.
    """

    path = Path(name=name, root=path)

    if force_download:
        os.system(f"rm -rf {path.get_main_folder()}")

    if not exists(path.get_data_json()):
        print("{}: Downloading dataset from HuggingFace Hub...".format(name))
        download_dataset(name)
        print(
            "{}: Download complete. File saved at {}".format(name, path.get_data_json())
        )

    return json.load(open(path.get_data_json(), "r"))


def download_dataset(name: str, root="data"):
    """Download a dataset from HuggingFace Hub. The name corresponds to the name of the dataset on HuggingFace Hub."""
    snapshot_download(
        repo_id="rouskinlab/" + name,
        repo_type="dataset",
        local_dir=Path(name, root).get_main_folder(),
        token=Env.get_hf_token(),
        allow_patterns=["data.json"],
    )


def name_from_path(datapath: str):
    if datapath.split("/")[-1] == "data.json":
        return datapath.split("/")[-2]
    return datapath.split("/")[-1].split(".")[0]


def clean_data(datapath: str):
    data = json.load(open(datapath, "r"))
    for ref, values in data.items():
        copy = values.copy()
        for k, v in values.items():
            if v is None or (isinstance(v, float) and np.isnan(v)):
                copy.pop(k)
        data[ref] = copy
    return data


def upload_dataset(
    datapath: str, exist_ok=False, commit_message: str = None, add_card=True, **kwargs
):
    api = HfApi()
    name = name_from_path(datapath)
    # data = clean_data(datapath)

    hf_token = Env.get_hf_token()

    api.create_repo(
        repo_id="rouskinlab/" + name,
        token=hf_token,
        exist_ok=exist_ok,
        private=True,
        repo_type="dataset",
    )

    api.upload_file(
        path_or_fileobj=datapath,
        path_in_repo="data.json",
        repo_id="rouskinlab/" + name,
        repo_type="dataset",
        token=hf_token,
        commit_message=commit_message,
        **kwargs,
    )

    if add_card:
        card = write_card(datapath)
        api.upload_file(
            path_or_fileobj=card,
            repo_id="rouskinlab/" + name,
            path_in_repo="README.md",
            repo_type="dataset",
            token=hf_token,
            commit_message=commit_message,
        )


def write_card(datapath):
    source = os.path.basename(datapath)
    date = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    name = name_from_path(datapath)

    out = """
---
license: mit
language:
  - en
tags:
  - chemistry
  - biology`
author: Silvi Rouskin
source: {}
date: {}
---


# Data types
""".format(
        source, date
    )

    data_type_count = {}
    data = json.load(open(datapath, "r"))
    for dp in data.values():
        for k, v in dp.items():
            if v is None or (isinstance(v, float) and np.isnan(v)):
                continue
            if k not in data_type_count:
                data_type_count[k] = 0
            data_type_count[k] += 1

    for k, v in data_type_count.items():
        out += f"""
- **{k}**: {v} datapoints"""

    # TODO add filtering report

    # dump
    path = Path(name=name, root=dirname(dirname(datapath)))
    os.makedirs(os.path.dirname(path.get_card()), exist_ok=True)
    with open(path.get_card(), "w") as f:
        f.write(out)

    if os.path.exists(path.get_conversion_report()):
        conversion = open(path.get_conversion_report(), "r").read()
        with open(path.get_card(), "a") as f:
            f.write("\n\n" + conversion)

    return path.get_card()
