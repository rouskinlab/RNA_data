
from .env import DATA_FOLDER
import os
from os.path import join

class PathDatafolder:
    """Path to files and folders of a datafolder of name `name`. The path to the data folder is `DATA_FOLDER`, which is defined in `env`.

    Parameters
    ----------

    name : str

    Returns
    -------

    PathDatafolder

    Example
    -------

    >>> path = PathDatafolder(name='my_test_datafolder_pytest')
    >>> path.name
    'my_test_datafolder_pytest'
    >>> print(path)
    PathDatafolder(name='my_test_datafolder_pytest')
    """

    def __init__(self, name, root = DATA_FOLDER) -> None:
        assert type(name) == str, f'name {name} is not a string'
        assert type(root) == str, f'root {root} is not a string'
        self.root = root
        self.name = name

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name='{self.name}')"

    def get_data_folder(self)->str:
        """Returns the path to the data folder."""
        return self.root

    def get_main_folder(self)->str:
        """Returns the path to the main folder."""
        return join(self.get_data_folder(), self.name)

    def get_dms_npy(self)->str:
        """Returns the path to the DMS npy file"""
        return join(self.get_main_folder(), 'dms.npy')

    def get_base_pairs_npy(self)->str:
        """Returns the path to the structure npy file"""
        return join(self.get_main_folder(), 'base_pairs.npy')

    def get_references_npy(self)->str:
        """Returns the path to the names npy file"""
        return join(self.get_main_folder(), 'references.npy')

    def get_sequences_npy(self)->str:
        """Returns the path to the sequences npy file"""
        return join(self.get_main_folder(), 'sequences.npy')

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

    def get_readme_file(self):
        """Returns the path to the readme file."""
        return join(self.get_main_folder(), "README.md")
