import os
from .util import source_env

HUGGINGFACE_TOKEN = os.environ['HUGGINGFACE_TOKEN']
DATA_FOLDER = os.environ['DATA_FOLDER']
DATA_FOLDER_TESTING = os.environ['DATA_FOLDER_TESTING']
RNASTRUCTURE_PATH = os.environ['RNASTRUCTURE_PATH']
RNASTRUCTURE_TEMP_FOLDER = os.environ['RNASTRUCTURE_TEMP_FOLDER']

for var in ["HUGGINGFACE_TOKEN", "DATA_FOLDER", "DATA_FOLDER_TESTING", "RNASTRUCTURE_PATH", "RNASTRUCTURE_TEMP_FOLDER"]:
    if len(os.environ[var]) == 0:
        raise ValueError(f"Environment variable {var} is not set. Please run `source env`.")

