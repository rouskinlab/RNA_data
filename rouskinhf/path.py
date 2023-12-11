from os.path import join
import os
from .env import Env


class Path:
    """Path to files and folders of a datafolder of name `name`. The path to the data folder is `DATA_FOLDER`, which is defined in `env`.

    Parameters
    ----------

    name : str

    Returns
    -------

    Path

    Example
    -------

    >>> path = Path(name='my_test_datafolder_pytest', root='data')
    >>> path.name
    'my_test_datafolder_pytest'
    >>> print(path)
    Path(name='my_test_datafolder_pytest')
    """

    def __init__(self, name, root) -> None:
        assert isinstance(name, str), f"name {name} is not a string"
        assert isinstance(root, str), f"root {root} is not a string"
        self.root = root
        self.name = name

    def make(self, force=False):
        """Creates the data folder."""
        if force:
            self.clear()
        os.system(f"mkdir -p {self.get_main_folder()}")

    def clear(self):
        """Clears the data folder."""
        os.system(f"rm -rf {self.get_main_folder()}")

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name='{self.name}')"

    def get_data_folder(self) -> str:
        """Returns the path to the data folder."""
        return self.root

    def get_main_folder(self) -> str:
        """Returns the path to the main folder."""
        return join(self.get_data_folder(), self.name)

    def get_data_json(self) -> str:
        """Returns the path to the data.json file."""
        return join(self.get_main_folder(), "data.json")

    def get_card(self) -> str:
        """Returns the path to the README.md file."""
        return join(self.get_main_folder(), "README.md")

    def get_conversion_report(self) -> str:
        """Returns the path to the conversion report file."""
        return join(self.get_main_folder(), "conversion_report.txt")
