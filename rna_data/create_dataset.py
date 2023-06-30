
from typing import Any

from rna_data.config import DATA_FOLDER
from .path_dataset import PathDataset
from .config import DATA_FOLDER
from .datapoints import ListofDatapoints, write_list_of_datapoints_to_json
from .info_file import infoFileWriter
from .write_npy import write_dms_npy_from_json, write_structure_npy_from_json
import os
import numpy as np
from huggingface_hub import HfApi
from .config import HUGGINGFACE_TOKEN
from huggingface_hub import snapshot_download

GENERATE_NPY = False
PREDICT_STRUCTURE = False
PREDICT_DMS = False
ROUSKINLAB = 'rouskinlab/'

class DataFolder:
    """Create a datafolder from a fasta file, a json file or a folder of ct files.

    Parameters
    ----------

    path_in : str
        path_in to the fasta file, the json file or the folder of ct files.

    path_out : str
        path_out to the folder where the dataset is created.

    name : str, optional
        Name of the dataset. If None, the name is the name of the file or folder.

    predict_structure : bool, optional
        If True, the structure is predicted. Default is True.

    predict_dms : bool, optional
        If True, the dms is predicted. Default is True.

    Returns
    -------

    dataset : Dataset

    """

    def from_fasta(path_in, path_out=DATA_FOLDER, name = None, predict_structure = PREDICT_STRUCTURE, predict_dms = PREDICT_DMS, generate_npy = GENERATE_NPY):
        return CreateDatasetFromFasta(path_in, path_out, name, predict_structure, predict_dms, generate_npy)

    def from_dreem_output(path_in, path_out=DATA_FOLDER, name = None, predict_structure = PREDICT_STRUCTURE, generate_npy = GENERATE_NPY):
        return CreateDatasetFromDreemOutput(path_in, path_out, name, predict_structure, generate_npy)

    def from_json(path_in, path_out=DATA_FOLDER, name = None, generate_npy = GENERATE_NPY):
        return CreateDatasetFromJSON(path_in, path_out, name, generate_npy)

    def from_ct_folder(path_in, path_out=DATA_FOLDER, name = None, predict_dms = PREDICT_DMS, generate_npy = GENERATE_NPY):
        return CreateDatasetFromCTfolder(path_in, path_out, name, predict_dms, generate_npy)

    def from_huggingface(name, path_out=DATA_FOLDER):
        return LoadDatasetFromHF(name, path_out)



class CreateDatasetTemplate(PathDataset):

    def __init__(self, path_in, path_out, name, source, predict_structure, predict_dms) -> None:
        super().__init__(name, path_out)
        self.path_in = path_in
        self.path_out = path_out
        self.datapoints = []
        self.api = HfApi(token=HUGGINGFACE_TOKEN)
        self.predict_structure = predict_structure
        self.predict_dms = predict_dms

        # Set name
        if name is None:
            name = path_in.replace('\\','/').split('/')[-1].split('.')[0]
        self.name = name

        # move path_in to source folder
        os.makedirs(self.get_source_folder(), exist_ok=True)
        os.system(f'cp -fr {path_in} {self.get_source_folder()}')

        # Write info file
        infoFileWriter(source=source, dataset=self).write()


    def __repr__(self) -> str:
        return f"{self.__class__.__name__} @{self.get_main_folder()}"

    def create_repo(self, exist_ok=False, private=True):
        self.api.create_repo(
            repo_id=ROUSKINLAB+self.name,
            token=HUGGINGFACE_TOKEN,
            exist_ok=exist_ok,
            private=private,
        )

    def upload_folder(self, revision = 'main', commit_message=None, commit_description=None, multi_commits=False, run_as_future=False, **kwargs):
        future = self.api.upload_folder(
            repo_id=ROUSKINLAB+self.name,
            folder_path=self.get_main_folder(),
            token=HUGGINGFACE_TOKEN,
            revision=revision,
            commit_message=commit_message,
            commit_description=commit_description,
            multi_commits=multi_commits,
            run_as_future=run_as_future,
            **kwargs
        )

        if run_as_future:
            return future

    def generate_npy(self):

        if self.predict_structure:
            write_structure_npy_from_json(self)

        if self.predict_dms:
            write_dms_npy_from_json(self)


class CreateDatasetFromJSON(CreateDatasetTemplate):

    def __init__(self, path_in, name, predict_structure, predict_dms, generate_npy) -> None:
        super().__init__(path_in, name)


class CreateDatasetFromDreemOutput(CreateDatasetTemplate):

    """Create a dataset from a dreem output file.

    Parameters
    ----------

    path_in : str
        path_in to the dreem output file.

    name : str, optional
        Name of the dataset. If None, the name is the name of the file or folder.

    predict_structure : bool, optional
        If True, the structure of the RNA is predicted. Default is True.

    generate_npy : bool, optional
        If True, the npy files are generated. Default is True.

    >>> dataset = DataFolder.from_dreem_output(path_in='data/dreem_output.json', generate_npy=True)
    >>> dataset.name
    'dreem_output'
    >>> print(dataset)
    CreateDatasetFromDreemOutput @data/datasets/dreem_output
    >>> os.path.isfile(dataset.get_json())
    True
    """

    def __init__(self, path_in, path_out, name, predict_structure, generate_npy) -> None:
        super().__init__(path_in, path_out, name, source = 'dreem_output', predict_structure = predict_structure, predict_dms = False)

        write_list_of_datapoints_to_json(
            path = self.get_json(),
            datapoints = ListofDatapoints.from_dreem_output(path_in, predict_structure = predict_structure)
        )

        if generate_npy:
            self.generate_npy()


class CreateDatasetFromFasta(CreateDatasetTemplate):

    """ Create a dataset from a fasta file.

    Parameters
    ----------

    path_in : str
        path_in to the dreem output file.

    name : str, optional
        Name of the dataset. If None, the name is the name of the file or folder.

    predict_structure : bool, optional
        If True, the structure of the RNA is predicted. Default is True.

    predict_dms : bool, optional
        If True, the dms of the RNA is predicted. Default is True.

    generate_npy : bool, optional
        If True, the npy files are generated. Default is True.

    >>> dataset = DataFolder.from_fasta(path_in='data/sequences.fasta', generate_npy=True)
    >>> dataset.name
    'sequences'
    >>> print(dataset)
    CreateDatasetFromFasta @data/datasets/sequences
    >>> os.path.isfile(dataset.get_json())
    True
    """

    def __init__(self, path_in, path_out, name, predict_structure, predict_dms, generate_npy) -> None:
        super().__init__(path_in, path_out, name, source = 'fasta', predict_structure = predict_structure, predict_dms = predict_dms)

        write_list_of_datapoints_to_json(
            path = self.get_json(),
            datapoints = ListofDatapoints.from_fasta(path_in, predict_structure = predict_structure, predict_dms = predict_dms)
        )

        if generate_npy:
            self.generate_npy()


class CreateDatasetFromCTfolder(CreateDatasetTemplate):

    """ Create a dataset from a folder of ct files.

    Parameters
    ----------

    path_in : str
        path_in to the dreem output file.

    name : str, optional
        Name of the dataset. If None, the name is the name of the file or folder.

    predict_dms : bool, optional
        If True, the dms of the RNA is predicted. Default is True.

    generate_npy : bool, optional
        If True, the npy files are generated. Default is True.

    >>> dataset = DataFolder.from_ct_folder(path_in='data/ct_files', generate_npy=True)
    >>> dataset.name
    'ct_files'
    >>> print(dataset)
    CreateDatasetFromCTfolder @data/datasets/ct_files
    >>> os.path.isfile(dataset.get_json())
    True
    """

    def __init__(self, path_in, path_out, name, predict_dms, generate_npy) -> None:
        super().__init__(path_in, path_out, name, source = 'ct', predict_structure = False, predict_dms = predict_dms)

        ct_files = [os.path.join(path_in, f) for f in os.listdir(path_in) if f.endswith('.ct')]
        write_list_of_datapoints_to_json(
            path = self.get_json(),
            datapoints = ListofDatapoints.from_ct(ct_files, predict_dms = predict_dms)
        )

        if generate_npy:
            self.generate_npy()


class LoadDatasetFromHF(PathDataset):
    """Load a dataset from HuggingFace.

    Parameters
    ----------

    name : str
        Name of the dataset.

    path_out : str
        Path to the folder where the dataset is saved.

    revision : str, optional
        Revision of the dataset. Default is 'main'.

    >>> dataset = LoadDatasetFromHF(name='dreem_output', path_out='data/datasets')
    >>> dataset.name
    'dreem_output'
    """

    def __init__(self, name, path_out, revision='main') -> None:
        super().__init__(name, path_out)

        # Download the dataset #TODO : check if the dataset is already downloaded
        snapshot_download(
            repo_id = ROUSKINLAB+self.name,
         #   repo_type='space',
            local_dir=self.get_main_folder(),
            revision=revision,
            token=HUGGINGFACE_TOKEN)

        assert os.path.isdir(self.get_main_folder()), f'No folder found in {self.get_main_folder()}'
        assert os.path.isfile(self.get_json()), f'No json file found in {self.get_main_folder()}'



