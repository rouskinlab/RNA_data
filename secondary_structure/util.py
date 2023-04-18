import pandas as pd
import os
import numpy as np

import sys

# Define the one-hot encodings for the sequences and structures
seq2int = {
        'A': 1,
        'C': 2,
        'G': 3,
        'T': 4,
        'U': 4,
        'Y': 5,
        'R': 6,
        'K': 7,
        'W': 8,
        'S': 9,
        'M': 10,
        'B': 11,
        'D': 12,
        'H': 13,
        'V': 14,
        'N': 15,
        'X': 0
    }

struct2int = {
    'f': 1,
    't': 2,
    'i': 3,
    'h': 4,
    'm': 5,
    's': 6,
    'X': 0
}

dot2int = {'.': 1, '(': 2, ')': 3, 'X': 0}
int2dot = ['X', '.', '(', ')']

def import_structure(path_to_structures=None, data_type = 'full', size=None, save=False, reload=True):
    """
    Import the secondary structure dataset and convert to integer encoding, with padding.

    Each row of the dataset contains a sequence and a structure.
    The sequences contains the nucleotides A, C, G, T, U, Y, R, K, W, S, M, B, D, H, V, N
    The structures is a dot-bracket notation e.g. (((.)))

    :param path_to_structures: Path to the secondary structure dataset in json format
    :param size: The number of datapoints to import
    :param save: Whether to save the dataset as a numpy array
    :param reload: Whether to reload the dataset from the numpy array
    :param test: Whether to import the test dataset

    :return: A tuple with the integer encoded sequences and structures in numpy arrays
    """

    assert data_type in ['full', 'test', 'train', 'binary'], "Type must be 'full', 'test', 'binary', or 'train'"

    # Paths to the dataset
    dirname = os.path.dirname(os.path.abspath(__file__))
    save_path = [os.path.join(dirname, 'dataset', data_type, 'processed_sequences.npy'),
                 os.path.join(dirname, 'dataset', data_type, 'processed_structures.npy')]

    if path_to_structures is None:
        path_to_structures = os.path.join(dirname, 'dataset', data_type, 'secondary_structure.json')

    # Import the dataset and check size
    df = pd.read_json(path_to_structures)
    if size is None:
        size = len(df)
    elif size > len(df):
        print("Requested size too large, using full dataset")
        size = len(df)

    # Check if the dataset is already saved
    if reload:
        if os.path.exists(save_path[0]) and os.path.exists(save_path[1]):
            print("Loading saved dataset")
            sequences = np.load(save_path[0])
            structures = np.load(save_path[1])

            if len(sequences) < size or len(structures) < size:
                print("Dataset too small, creating new one")
                return import_structure(path_to_structures, data_type=data_type, size=size, save=save, reload=False)
            else:
                idx = np.random.choice(len(sequences), size=size, replace=False)
                if save:
                    np.save(save_path[0], sequences[idx])
                    np.save(save_path[1], structures[idx])
                
                return sequences[idx], structures[idx]
        else:
            print("Dataset not found, creating new one")
            return import_structure(path_to_structures, data_type=data_type, size=size, save=save, reload=False)
        

    else:        
        
        # Get a random sample of the dataset
        idx = np.random.choice(len(df), size=size, replace=False)
        df = df.iloc[idx].reset_index(drop=True)

        # Get max sequence length
        max_seq_len = df['sequence'].str.len().max()

        # Init numpy arrays for sequences and structures
        sequences = np.zeros((len(df), max_seq_len))
        structures = np.zeros((len(df), max_seq_len))

        # Iterate over the rows of the dataframe
        for i, row in df.iterrows():

            # Print progress
            if i%1000:
                sys.stdout.write("Processing dataset: %d%%   \r" % (100*i/len(df)) )
                sys.stdout.flush()
            
            # Apply zero padding to each row (sequence and structure)
            row['sequence'] = row['sequence'].upper().ljust(max_seq_len, 'X')
            row['structure'] = row['structure'].ljust(max_seq_len, 'X')
            # row['structure'] = row['structure'] + [0]*(max_seq_len-len(row['structure']))

            # Integer encoding of the sequence, with padding as 'X' 
            for j, base in enumerate(row['sequence']):
                sequences[i,j] = seq2int[base]

            # One hot encoding of the structure, with padding as 'X' 
            for j, struct in enumerate(row['structure']):
                structures[i,j] = dot2int[struct]
                # structures[i,j] = struct

        # Save sequences and structures as numpy arrays
        if save:
            np.save(save_path[0], sequences)
            np.save(save_path[1], structures)

        
        return sequences, structures



if __name__ == '__main__':

    sequences, structures = import_structure(data_type='binary', save=True, reload=False)
    print("Loaded full dataset with shape: \n", sequences.shape, "\n", structures.shape)