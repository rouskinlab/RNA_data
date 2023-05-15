import os
import numpy as np
import sys
import pandas as pd
import json
from tqdm import tqdm

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'RNAstructure'))
from src.util import fastaToDict
from RNAstructure import predictFromSequence


def generate_pairing_matrix(path_to_fasta, path_to_RNAStructure=''):

    # Read the fasta file and get the list of sequences
    data = fastaToDict(path_to_fasta)
    ref_todelete = []
    for ref in tqdm(data.keys()):

        # Check if the sequence contains only A, C, G, U/T, discard it from the dataset otherwise
        bases = ['A','U','C','G','T','N','a','u','c','g','t','n']
        if not all(c in bases for c in data[ref]['sequence']):
            print(f'WARNING: sequence {ref} contains non-canonical bases, skipping it')
            ref_todelete.append(ref)
            continue

        prediction = predictFromSequence(data[ref]['sequence'], predict_pairs=True, predict_structure=False, predict_pairing_probability=False, rnastructure_path=path_to_RNAStructure)

        paired_bases = np.array(prediction['list_pairs'])-1
        data[ref]['paired_bases'] = paired_bases.tolist()

        # Check constraints -> !! Could be moved to dataset analysis !!
        pairing_matrix = np.zeros((len(data[ref]['sequence']), len(data[ref]['sequence'])))
        if len(paired_bases) > 0:
            pairing_matrix[paired_bases[:,0], paired_bases[:,1]] = 1
            pairing_matrix[paired_bases[:,1], paired_bases[:,0]] = 1

        assert np.all(pairing_matrix == pairing_matrix.T), 'ERROR: pairing matrix is not symmetric'
        assert ((np.sum(pairing_matrix, axis=1)==1) | (np.sum(pairing_matrix, axis=1)==0)).all(), 'ERROR: pairing matrix has more than one pairing per base'


    # Delete sequences with non-canonical bases
    for ref in ref_todelete:
        del data[ref]
        
    return data


if __name__ == '__main__':

    dir_name = os.path.dirname(os.path.abspath(__file__))

    # Get path to fasta file
    if len(sys.argv) >= 2:
        path_to_fasta = sys.argv[1]
    else:
        path_to_fasta = os.path.join(dir_name, '..', 'sequence_dataset', 'sequences_half.fasta')

    # Get the secondary structure sequences as pairing matrix
    data = generate_pairing_matrix(path_to_fasta, path_to_RNAStructure='/Users/alberic/RNAstructure/exe/')

    # Save the dataframes as json
    save_dir = os.path.join(dir_name, 'dataset', 'pairing_matrix')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    # Write to json file
    with open(os.path.join(save_dir, 'secondary_structure.json'), "w") as outfile:
        json.dump(data, outfile, indent=2)