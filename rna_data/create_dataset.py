
from typing import Any
from .path_dataset import PathDataset
from .config import DATA_FOLDER
from .datapoints import ListofDatapoints, write_list_of_datapoints_to_json
from .info_file import infoFileWriter
from .write_npy import write_npy_from_json
import os
import numpy as np

GENERATE_NPY = False
PREDICT_STRUCTURE = False
PREDICT_DMS = False

class CreateDataset:
    """Create a dataset from a fasta file, a json file or a folder of ct files.

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


class CreateDatasetTemplate(PathDataset):

    def __init__(self, path_in, path_out, name, source) -> None:
        super().__init__(name, path_out)
        self.path_in = path_in
        self.path_out = path_out
        self.datapoints = []

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

    def push_to_hub(self, overwrite=False):
        pass

    def generate_npy(self):
        write_npy_from_json(
            json_path=self.get_json(),
            npy_path=self.get_npy_folder()
        )


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

    >>> dataset = CreateDataset.from_dreem_output(path_in='data/dreem_output.json', generate_npy=True)
    >>> dataset.name
    'dreem_output'
    >>> print(dataset)
    CreateDatasetFromDreemOutput @data/datasets/dreem_output
    >>> os.listdir(dataset.get_npy_folder())
    ['placeholder.npy']
    >>> os.path.isfile(dataset.get_json())
    True
    """

    def __init__(self, path_in, path_out, name, predict_structure, generate_npy) -> None:
        super().__init__(path_in, path_out, name, source = 'dreem_output')

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

    >>> dataset = CreateDataset.from_fasta(path_in='data/sequences.fasta', generate_npy=True)
    >>> dataset.name
    'sequences'
    >>> print(dataset)
    CreateDatasetFromFasta @data/datasets/sequences
    >>> os.listdir(dataset.get_npy_folder())
    ['placeholder.npy']
    >>> os.path.isfile(dataset.get_json())
    True
    """

    def __init__(self, path_in, path_out, name, predict_structure, predict_dms, generate_npy) -> None:
        super().__init__(path_in, path_out, name, source = 'fasta')

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

    >>> dataset = CreateDataset.from_ct_folder(path_in='data/ct_files', generate_npy=True)
    >>> dataset.name
    'ct_files'
    >>> print(dataset)
    CreateDatasetFromCTfolder @data/datasets/ct_files
    >>> os.listdir(dataset.get_npy_folder())
    ['placeholder.npy']
    >>> os.path.isfile(dataset.get_json())
    True
    """

    def __init__(self, path_in, path_out, name, predict_dms, generate_npy) -> None:
        super().__init__(path_in, path_out, name, source = 'ct')

        ct_files = [os.path.join(path_in, f) for f in os.listdir(path_in) if f.endswith('.ct')]
        write_list_of_datapoints_to_json(
            path = self.get_json(),
            datapoints = ListofDatapoints.from_ct(ct_files, predict_dms = predict_dms)
        )

        if generate_npy:
            self.generate_npy()

