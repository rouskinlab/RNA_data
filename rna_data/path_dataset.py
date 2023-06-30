
from .config import DATA_FOLDER
import os


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
    >>> path_dataset.get_npy_folder()
    'data/datasets/my_test_dataset_pytest/npy'
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
        return f"{self.root}/{self.name}"

    def get_npy_folder(self)->str:
        """Returns the path to the npy folder."""
        return f"{self.root}/{self.name}/npy"

    def get_npy_files(self)->list:
        """Returns the list of npy files, if the npy folder exists."""
        return [os.path.join(self.get_npy_folder(), f) for f in os.listdir(self.get_npy_folder())]

    def get_json(self)->str:
        """Returns the path to the json."""
        return f"{self.root}/{self.name}/data.json"

    def get_source_folder(self):
        """Returns the path to the source folder."""
        return f"{self.root}/{self.name}/source"

    def get_source_files(self):
        """Returns the list of source files, if the source folder exists."""
        return [os.path.join(self.get_source_folder(), f) for f in os.listdir(self.get_source_folder())]

    def get_info_file(self):
        """Returns the path to the info file."""
        return f"{self.root}/{self.name}/info.json"
