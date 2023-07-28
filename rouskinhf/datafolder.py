
from typing import Any

import json
from .env import DATA_FOLDER

from .path import PathDatafolder
from .env import DATA_FOLDER, HUGGINGFACE_TOKEN
from .list_datapoints import ListofDatapoints
from .info_file import infoFileWriter
import os
import numpy as np
from huggingface_hub import HfApi
from huggingface_hub import snapshot_download


GENERATE_NPY = True
PREDICT_STRUCTURE = False
PREDICT_DMS = False
ROUSKINLAB = 'rouskinlab/'

class DataFolderTemplate(PathDatafolder):

    def __init__(self, name, path_out) -> None:
        super().__init__(name, path_out)
        self.name = name
        self.datapoints = ListofDatapoints([], verbose=False)

    def generate_npy(self):
        """Generate the npy files:
            - reference.npy
            - sequence.npy
            - base_pairs.npy
            - dms.npy

        Each file corresponds to a key in the datapoints dictionary contained in the json file.
        For example, if the 1st line of the json file is:
        {"RF02271.fa.csv_1":{"sequence": "AAACCCCAGUAGGGGUUU", "paired_bases": [[6, 11], [5, 12], [4, 13], [3, 14], [2, 15], [1, 16], [0, 17]], "dms": [0.3317625732102161, 0.8619477949493498, 0.9849505502538177, 0.9997828311396947, 0.9999562120800795, 0.9999454042235505, 0.9970246205212348, 3.298876877182004e-10, 1.0113344525791075e-08, 9.30187856307975e-09, 2.1481183689942838e-10, 0.996953491189411, 0.9999737756512344, 0.9999846271614127, 0.9997971647502685, 0.9848757245213886, 0.8618760196213241, 0.3319091664146036]}}
        then the 1st line of the reference.npy file will be:
        AAACCCCAGUAGGGGUUU
        """

        if not len(self.datapoints.datapoints):
            return None

        self.datapoints.to_reference_npy(self.get_references_npy())
        self.datapoints.to_sequence_npy(self.get_sequences_npy())

        d1 = self.datapoints.datapoints[0]

        if hasattr(d1, 'paired_bases'):
            self.datapoints.to_base_pairs_npy(self.get_base_pairs_npy())

        if hasattr(d1, 'dms'):
            self.datapoints.to_dms_npy(self.get_dms_npy())

class CreateDatafolderTemplate(DataFolderTemplate):

    def __init__(self, path_in, path_out, name, source, predict_structure, predict_dms) -> None:
        name = self._set_name(name, path_in)
        super().__init__(name, path_out)
        self.name = name
        self.path_in = path_in
        self.path_out = path_out
        self.api = HfApi(token=HUGGINGFACE_TOKEN)
        self.predict_structure = predict_structure
        self.predict_dms = predict_dms

        # move path_in to source folder
        os.makedirs(self.get_source_folder(), exist_ok=True)
        os.system(f'cp -fr {path_in} {self.get_source_folder()}')

        # Write info file
        self.infofile = infoFileWriter(source=source, datafolder=self, predict_dms=predict_dms, predict_structure=predict_dms)


    def __repr__(self) -> str:
        return f"{self.__class__.__name__} @{self.get_main_folder()}"

    def _set_name(self, name, path_in):
        if name is None:
            name = path_in.replace('\\','/').split('/')[-1].split('.')[0]
        return name

    def create_repo(self, exist_ok=False, private=True):

        """Create a repo on huggingface.co.

        Parameters
        ----------

        exist_ok : bool, optional
            If True, the repo is created even if it already exists. Default is False.

        private : bool, optional
            If True, the repo is private. Default is True.

        Examples
        --------

        >>> datafolder = DataFolder.from_dreem_output(path_in='data/input_files_for_testing/dreem_output.json', verbose=False)
        >>> datafolder.create_repo(exist_ok=True)
        """

        self.api.create_repo(
            repo_id=ROUSKINLAB+self.name,
            token=HUGGINGFACE_TOKEN,
            exist_ok=exist_ok,
            private=private,
            repo_type="dataset",
        )

    def upload_folder(self, revision = 'main', commit_message=None, commit_description=None, multi_commits=False, run_as_future=False, **kwargs):

        """Upload the datafolder to huggingface.co. Only the json file, the README.md and the source files are uploaded.

        Parameters
        ----------

        revision : str, optional
            The revision of the repo (i.e. the branch on the HuggingFace hub). Default is 'main'.

        commit_message : str, optional
            The commit message. Default is None.

        commit_description : str, optional
            The commit description. Default is None.

        multi_commits : bool, optional
            Useful for large files. If True, the upload is done in multiple commits. Default is False.

        run_as_future : bool, optional
            If True, the upload is run asynchonously. The state of the upload can be checked with the returned future. See sample code. Default is False.

        Examples
        --------

        >>> datafolder = DataFolder.from_dreem_output(path_in='data/input_files_for_testing/dreem_output.json', verbose=False)
        >>> datafolder.create_repo(exist_ok=True)
        >>> datafolder.upload_folder('main', 'pytest commit', 'A commit generated by pytest for the upload_folder method')
        >>> future = datafolder.upload_folder(run_as_future=True)
        >>> future.result() # wait for the upload to finish
        'https://huggingface.co/datasets/rouskinlab/dreem_output/tree/main/'
        >>> future.done() # check if the upload is done
        True
        """

        future = self.api.upload_folder(
            repo_id=ROUSKINLAB+self.name,
            folder_path=self.get_main_folder(),
            repo_type="dataset",
            token=HUGGINGFACE_TOKEN,
            revision=revision,
            commit_message=commit_message,
            commit_description=commit_description,
            multi_commits=multi_commits,
            run_as_future=run_as_future,
            allow_patterns=["source/*", "info.json", "data.json", "README.md"],
            **kwargs
        )

        if run_as_future:
            return future


    def dump_datapoints(self, generate_npy):
        """Dump the datapoints to a json file."""

        self.datapoints.to_json(self.get_json())

        if generate_npy:
            self.generate_npy()



class CreateDatafolderFromDreemOutput(CreateDatafolderTemplate):

    """Create a datafolder from a dreem output file.

    Parameters
    ----------

    path_in : str
        path_in to the dreem output file.

    name : str, optional
        Name of the datafolder. If None, the name is the name of the file or folder.

    predict_structure : bool, optional
        If True, the structure of the RNA is predicted. Default is True.

    generate_npy : bool, optional
        If True, the npy files are generated. Default is True.

    tqdm : bool, optional
        If True, a progress bar is displayed. Default is True.


    Examples
    --------

    >>> datafolder = DataFolder.from_dreem_output(path_in='data/input_files_for_testing/dreem_output.json', generate_npy=True)
    Over a total of 13 datapoints, there are:
        - 12 valid datapoints
        - 1 invalid datapoints (ex: sequence with non-regular characters)
        - 0 datapoints with the same reference
        - 0 duplicate sequences with the same structure / dms
        - 0 duplicate sequences with different structure / dms
    >>> datafolder.name
    'dreem_output'
    >>> print(str(datafolder).split(' ')[0])
    CreateDatafolderFromDreemOutput
    >>> os.path.isfile(datafolder.get_json())
    True
    """

    def __init__(self, path_in, path_out, name, predict_structure, generate_npy, tqdm=True, verbose=True) -> None:
        super().__init__(path_in, path_out, name, source = 'dreem_output', predict_structure = predict_structure, predict_dms = False)

        self.datapoints = ListofDatapoints.from_dreem_output(path_in, predict_structure = predict_structure, tqdm=tqdm, verbose=verbose)
        self.dump_datapoints(generate_npy)
        self.infofile.add_filtering_report(self.datapoints.filtering_report).write()

class CreateDatafolderFromFasta(CreateDatafolderTemplate):

    """ Create a datafolder from a fasta file.

    Parameters
    ----------

    path_in : str
        path_in to the dreem output file.

    name : str, optional
        Name of the datafolder. If None, the name is the name of the file or folder.

    predict_structure : bool, optional
        If True, the structure of the RNA is predicted. Default is True.

    predict_dms : bool, optional
        If True, the dms of the RNA is predicted. Default is True.

    generate_npy : bool, optional
        If True, the npy files are generated. Default is True.


    Examples
    --------

    >>> datafolder = DataFolder.from_fasta(path_in='data/input_files_for_testing/sequences.fasta', generate_npy=True)
    Over a total of 5 datapoints, there are:
        - 2 valid datapoints
        - 1 invalid datapoints (ex: sequence with non-regular characters)
        - 1 datapoints with the same reference
        - 1 duplicate sequences
    >>> datafolder.name
    'sequences'
    >>> print(str(datafolder).split(' ')[0])
    CreateDatafolderFromFasta
    >>> os.path.isfile(datafolder.get_json())
    True
    """

    def __init__(self, path_in, path_out, name, predict_structure, predict_dms, generate_npy, tqdm, verbose=True) -> None:
        super().__init__(path_in, path_out, name, source = 'fasta', predict_structure = predict_structure, predict_dms = predict_dms)

        self.datapoints = ListofDatapoints.from_fasta(path_in, predict_structure = predict_structure, predict_dms = predict_dms, tqdm=tqdm, verbose=verbose)
        self.dump_datapoints(generate_npy)
        self.infofile.add_filtering_report(self.datapoints.filtering_report).write()



class CreateDatafolderFromCTfolder(CreateDatafolderTemplate):

    """ Create a datafolder from a folder of ct files.

    Parameters
    ----------

    path_in : str
        path_in to the dreem output file.

    name : str, optional
        Name of the datafolder. If None, the name is the name of the file or folder.

    predict_dms : bool, optional
        If True, the dms of the RNA is predicted. Default is True.

    generate_npy : bool, optional
        If True, the npy files are generated. Default is True.

    tqdm : bool, optional
        If True, a progress bar is displayed. Default is True.

    Examples
    --------

    >>> datafolder = DataFolder.from_ct_folder(path_in='data/input_files_for_testing/ct_files', generate_npy=True)
    Over a total of 5 datapoints, there are:
        - 1 valid datapoints
        - 1 invalid datapoints (ex: sequence with non-regular characters)
        - 0 datapoints with the same reference
        - 1 duplicate sequences with the same structure / dms
        - 2 duplicate sequences with different structure / dms
    >>> datafolder.name
    'ct_files'
    >>> print(str(datafolder).split(' ')[0])
    CreateDatafolderFromCTfolder
    >>> os.path.isfile(datafolder.get_json())
    True
    """

    def __init__(self, path_in, path_out, name, predict_dms, generate_npy, tqdm=True, verbose=True) -> None:
        super().__init__(path_in, path_out, name, source = 'ct', predict_structure = False, predict_dms = predict_dms)

        ct_files = [os.path.join(path_in, f) for f in os.listdir(path_in) if f.endswith('.ct') or f.endswith('.txt')]
        self.datapoints = ListofDatapoints.from_ct(ct_files, predict_dms = predict_dms, tqdm=tqdm, verbose=verbose)
        self.dump_datapoints(generate_npy)
        self.infofile.add_filtering_report(self.datapoints.filtering_report).write()


class LoadDatafolder(DataFolderTemplate):

    def __init__(self, name, root) -> None:
        super().__init__(name, root)


class LoadDatafolderFromHF(LoadDatafolder):
    """Load a datafolder from HuggingFace.

    Parameters
    ----------

    name : str
        Name of the datafolder.

    path_out : str
        Path to the folder where the datafolder is saved.

    revision : str, optional
        Revision of the datafolder. Default is 'main'.


    Examples
    --------

    >>> datafolder = LoadDatafolderFromHF(name='for_testing', path_out='data/datafolders')
    Over a total of 2 datapoints, there are:
        - 2 valid datapoints
        - 0 invalid datapoints (ex: sequence with non-regular characters)
        - 0 datapoints with the same reference
        - 0 duplicate sequences with the same structure / dms
        - 0 duplicate sequences with different structure / dms
    >>> datafolder.name
    'for_testing'
    """

    def __init__(self, name, path_out, revision='main', tqdm=True, verbose=True) -> None:
        super().__init__(name, path_out)

        # Download the datafolder #TODO : check if the datafolder is already downloaded
        snapshot_download(
            repo_id = ROUSKINLAB+self.name,
            repo_type='dataset',
            local_dir=self.get_main_folder(),
            revision=revision,
            token=HUGGINGFACE_TOKEN,
            allow_patterns=['info.json', 'data.json']
            )

        assert os.path.isdir(self.get_main_folder()), f'No folder found in {self.get_main_folder()}'
        assert os.path.isfile(self.get_json()), f'No json file found in {self.get_main_folder()}'

        self.datapoints = ListofDatapoints.from_json(self.get_json(), tqdm=tqdm, verbose=verbose)


class LoadDatafolderFromLocal(LoadDatafolder):
    """Load a datafolder from local.

    Parameters
    ----------

    name : str
        Name of the datafolder.

    path : str
        Path to the folder where the datafolder is saved. Defaut is set in the env variable DATAFOLDER_PATH.

    Examples
    --------

    >>> datafolder = LoadDatafolderFromLocal(name='for_testing', path='data/datafolders', generate_npy=True)
    Over a total of 2 datapoints, there are:
        - 2 valid datapoints
        - 0 invalid datapoints (ex: sequence with non-regular characters)
        - 0 datapoints with the same reference
        - 0 duplicate sequences with the same structure / dms
        - 0 duplicate sequences with different structure / dms
    >>> datafolder.name
    'for_testing'
    """

    def __init__(self, name, path, generate_npy, tqdm=True, verbose=True) -> None:
        super().__init__(name, path)

        assert os.path.isdir(self.get_main_folder()), f'No folder found in {self.get_main_folder()}'
        assert os.path.isfile(self.get_json()), f'No json file found in {self.get_main_folder()}'

        if generate_npy:
            self.datapoints = ListofDatapoints.from_json(self.get_json(), tqdm=tqdm, verbose=verbose)
            self.generate_npy()

class DataFolder:
    """Create a datafolder from a fasta file, a json file or a folder of ct files.

    Parameters
    ----------

    path_in : str
        path_in to the fasta file, the json file or the folder of ct files.

    path_out : str
        path_out to the folder where the datafolder is created.

    name : str, optional
        Name of the datafolder. If None, the name is the name of the file or folder.

    predict_structure : bool, optional
        If True, the structure is predicted. Default is True.

    predict_dms : bool, optional
        If True, the dms is predicted. Default is True.

    tqdm : bool, optional
        If True, a progress bar is displayed. Default is True.

    verbose : bool, optional
        If True, the filtering report is displayed. Default is True.

    """

    def from_fasta(path_in, path_out=DATA_FOLDER, name = None, predict_structure = PREDICT_STRUCTURE, predict_dms = PREDICT_DMS, generate_npy = GENERATE_NPY, tqdm=True, verbose=True)->CreateDatafolderFromFasta:
        """Create a datafolder from a fasta file. See CreateDatafolderFromFasta for more details."""
        return CreateDatafolderFromFasta(path_in, path_out, name, predict_structure, predict_dms, generate_npy, tqdm=tqdm, verbose=verbose)

    def from_dreem_output(path_in, path_out=DATA_FOLDER, name = None, predict_structure = PREDICT_STRUCTURE, generate_npy = GENERATE_NPY, tqdm=True, verbose=True)->CreateDatafolderFromDreemOutput:
        """Create a datafolder from a dreem output file. See CreateDatafolderFromDreemOutput for more details."""
        return CreateDatafolderFromDreemOutput(path_in, path_out, name, predict_structure, generate_npy, tqdm=tqdm, verbose= verbose)

    def from_local(name, path=DATA_FOLDER, tqdm=True, verbose=True, generate_npy=False)->LoadDatafolderFromLocal:
        """Load a datafolder from local. See LoadDatafolderFromLocal for more details."""
        return LoadDatafolderFromLocal(name, path, tqdm=tqdm, verbose=verbose, generate_npy=generate_npy)

    def from_ct_folder(path_in, path_out=DATA_FOLDER, name = None, predict_dms = PREDICT_DMS, generate_npy = GENERATE_NPY, tqdm=True, verbose=True)->CreateDatafolderFromCTfolder:
        """Create a datafolder from a folder of ct files. See CreateDatafolderFromCTfolder for more details."""
        return CreateDatafolderFromCTfolder(path_in, path_out, name, predict_dms, generate_npy, tqdm=tqdm, verbose= verbose)

    def from_huggingface(name, path_out=DATA_FOLDER, tqdm=True, verbose=True)->LoadDatafolderFromHF:
        """Load a datafolder from HuggingFace. See LoadDatafolderFromHF for more details."""
        return LoadDatafolderFromHF(name, path_out, tqdm=tqdm, verbose=verbose)
