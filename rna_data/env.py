import os
from .util import source_env

path_env = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'env')

try:
    vars_env = source_env(path_env)
except FileNotFoundError:
    print("No env file found. Using environment variables.")
    vars_env = {}

HUGGINGFACE_TOKEN = os.environ.get('HUGGINGFACE_TOKEN')
DATA_FOLDER = os.environ.get('DATA_FOLDER')
DATA_FOLDER_TESTING = os.environ.get('DATA_FOLDER_TESTING')
RNASTRUCTURE_PATH = os.environ.get('RNASTRUCTURE_PATH')
RNASTRUCTURE_TEMP_FOLDER = os.environ.get('RNASTRUCTURE_TEMP_FOLDER')

for var in ["HUGGINGFACE_TOKEN", "DATA_FOLDER", "DATA_FOLDER_TESTING", "RNASTRUCTURE_PATH", "RNASTRUCTURE_TEMP_FOLDER"]:
    if exec(var) == None and var in vars_env and vars_env[var] != None:
        exec(f"{var} = vars_env['{var}']") # read from env file directly
    else:
        raise ValueError(f"Environment variable {var} is not set. Please run `source env`.")

