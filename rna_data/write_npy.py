
import os

def write_npy_from_json(json_path, npy_path, overwrite=False):
    assert os.path.isfile(json_path), f'json_path {json_path} does not exist'
    os.makedirs(npy_path, exist_ok=True)
    # TODO
    with open(os.path.join(npy_path, 'placeholder.npy'), 'w') as f:
        f.write('TODO')
