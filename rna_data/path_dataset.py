
from .config import DATA_FOLDER
import os
from os.path import join

class PathDataset:
    """Path to files and folders of a dataset of name `name`. The path to the data folder is `DATA_FOLDER`, which is defined in `config.py`.

    Parameters
    ----------

    name : str

    Returns
    -------

    path_dataset : PathDataset

    Example
    -------

    >>> path_dataset = PathDataset(name='my_test_dataset_pytest')
    >>> path_dataset.name
    'my_test_dataset_pytest'
    >>> print(path_dataset)
    PathDataset(name='my_test_dataset_pytest')
    >>> path_dataset.get_structure_npy()
    'data/datasets/my_test_dataset_pytest/structure.npy'
    """

    def __init__(self, name, root = DATA_FOLDER) -> None:
        self.root = root
        self.name = name

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name='{self.name}')"

    def get_data_folder(self)->str:
        """Returns the path to the data folder."""
        return self.root

    def get_main_folder(self)->str:
        """Returns the path to the main folder."""
        return join(self.get_data_folder(),self.name)

    def get_dms_npy(self)->str:
        """Returns the path to the DMS npy file"""
        return join(self.get_main_folder(), 'dms.npy')

    def get_structure_npy(self)->str:
        """Returns the path to the structure npy file"""
        return join(self.get_main_folder(), 'structure.npy')

    def get_json(self)->str:
        """Returns the path to the json."""
        return join(self.get_main_folder(), "data.json")

    def get_source_folder(self):
        """Returns the path to the source folder."""
        return join(self.get_main_folder(), "source")

    def get_source_files(self):
        """Returns the list of source files, if the source folder exists."""
        return [os.path.join(self.get_source_folder(), f) for f in os.listdir(self.get_source_folder())]

    def get_info_file(self):
        """Returns the path to the info file."""
        return join(self.get_main_folder(), "info.json")
