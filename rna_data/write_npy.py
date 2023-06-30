
import os

def write_structure_npy_from_json(dataset):
    assert os.path.isfile(dataset.get_json()), f'json_path {dataset.get_json()} does not exist'
    # TODO
    with open(dataset.get_structure_npy(), 'w') as f:
        f.write('TODO')

def write_dms_npy_from_json(dataset):
    assert os.path.isfile(dataset.get_json()), f'json_path {dataset.get_json()} does not exist'
    # TODO
    with open(dataset.get_dms_npy(), 'w') as f:
        f.write('TODO')
