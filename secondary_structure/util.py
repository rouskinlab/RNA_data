import pandas as pd
import os
import numpy as np

import sys

# Define the one-hot encodings for the sequences and structures
seq2vec = {
        'A': [1, 0, 0, 0, 0],
        'C': [0, 1, 0, 0, 0],
        'G': [0, 0, 1, 0, 0],
        'T': [0, 0, 0, 1, 0],
        'U': [0, 0, 0, 1, 0],
        'Y': [0, 1, 0, 1, 0],
        'R': [1, 0, 1, 0, 0],
        'K': [0, 0, 1, 1, 0],
        'W': [1, 0, 0, 1, 0],
        'S': [0, 1, 1, 0, 0],
        'M': [1, 0, 1, 0, 0],
        'B': [0, 1, 1, 1, 0],
        'D': [1, 0, 1, 1, 0],
        'H': [1, 1, 0, 1, 0],
        'V': [1, 1, 1, 0, 0],
        'N': [1, 1, 1, 1, 0],
        'X': [0, 0, 0, 0, 1]
    }

struct2vec = {
    'f': [1, 0, 0, 0, 0, 0, 0],
    't': [0, 1, 0, 0, 0, 0, 0],
    'i': [0, 0, 1, 0, 0, 0, 0],
    'h': [0, 0, 0, 1, 0, 0, 0],
    'm': [0, 0, 0, 0, 1, 0, 0],
    's': [0, 0, 0, 0, 0, 1, 0],
    'X': [0, 0, 0, 0, 0, 0, 1]
}

def import_structure(path_to_structures=None, type = 'full', size=None, save=False, reload=True):
    """
    Import the secondary structure dataset and convert to one-hot encoding.

    Each row of the dataset contains a sequence and a structure.
    The sequences contains A, U, C, G and Y, R, N bases.
    The structures contains f, t, i, h, m, or s characters.

    :param path_to_structures: Path to the secondary structure dataset in json format
    :param size: The number of datapoints to import
    :param save: Whether to save the dataset as a numpy array
    :param reload: Whether to reload the dataset from the numpy array
    :param test: Whether to import the test dataset

    :return: A tuple with the one-hot encoded sequences and structures in numpy arrays
    """

    assert type in ['full', 'test', 'train'], "Type must be 'full', 'test' or 'train'"

    # Paths to the dataset
    dirname = os.path.dirname(os.path.abspath(__file__))
    save_path = [os.path.join(dirname, 'dataset', type, 'processed_sequences.npy'),
                 os.path.join(dirname, 'dataset', type, 'processed_structures.npy')]

    if path_to_structures is None:
        path_to_structures = os.path.join(dirname, 'dataset', type, 'secondary_structure.json')

    # Import the dataset and checck size
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
                return import_structure(path_to_structures, type=type, size=size, save=save, reload=False)
            else:
                if save:
                    np.save(save_path[0], sequences[:size])
                    np.save(save_path[1], structures[:size])
                
                idx = np.random.choice(len(sequences), size=size, replace=False)
                return sequences[idx], structures[idx]
        else:
            print("Dataset not found, creating new one")
            return import_structure(path_to_structures, type=type, size=size, save=save, reload=False)
        

    else:        
        
        # Get a random sample of the dataset
        idx = np.random.choice(len(df), size=size, replace=False)
        df = df.iloc[idx].reset_index(drop=True)

        # Get max sequence length
        max_seq_len = df['sequence'].str.len().max()

        # Init numpy arrays for sequences and structures
        sequences = np.zeros((len(df), len(seq2vec['A']), max_seq_len))
        structures = np.zeros((len(df), len(struct2vec['f']), max_seq_len))

        # Iterate over the rows of the dataframe
        for i, row in df.iterrows():

            # Print progress
            if i%1000:
                sys.stdout.write("Processing dataset: %d%%   \r" % (100*i/len(df)) )
                sys.stdout.flush()
            
            # Apply zero padding to each row (sequence and structure)
            row['sequence'] = row['sequence'].upper().ljust(max_seq_len, 'X')
            row['structure'] = row['structure'].lower().ljust(max_seq_len, 'X') 

            # One hot encoding of the sequence, with padding as 'X' 
            for j, base in enumerate(row['sequence']):
                sequences[i, :, j] = seq2vec[base]

            # One hot encoding of the structure, with padding as 'X' 
            for j, struct in enumerate(row['structure']):
                structures[i, :, j] = struct2vec[struct]

        # Save sequences and structures as numpy arrays
        if save:
            np.save(save_path[0], sequences)
            np.save(save_path[1], structures)

        
        return sequences, structures



if __name__ == '__main__':

    sequences, structures = import_structure(type='full', save=True, reload=False)
    print("Loaded full dataset with shape: \n", sequences.shape, "\n", structures.shape)

    sequences, structures = import_structure(type='train', save=True, reload=False)
    print("Loaded train dataset with shape: \n", sequences.shape, "\n", structures.shape)

    sequences, structures = import_structure(type='test', save=True, reload=False)
    print("Loaded test dataset with shape: \n", sequences.shape, "\n", structures.shape)
