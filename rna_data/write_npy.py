
import os

def write_structure_npy_from_json(datafolder):
    assert os.path.isfile(datafolder.get_json()), f'json_path {datafolder.get_json()} does not exist'
    # TODO
    with open(datafolder.get_structure_npy(), 'w') as f:
        f.write('TODO')

def write_dms_npy_from_json(datafolder):
    assert os.path.isfile(datafolder.get_json()), f'json_path {datafolder.get_json()} does not exist'
    # TODO
    with open(datafolder.get_dms_npy(), 'w') as f:
        f.write('TODO')
