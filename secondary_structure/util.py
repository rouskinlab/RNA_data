import pandas as pd
import os
import numpy as np

import sys

# Define the one-hot encodings for the sequences and structures
seq2int = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3,
        'U': 3,
        'Y': 4,
        'R': 5,
        'K': 6,
        'W': 7,
        'S': 8,
        'M': 9,
        'B': 10,
        'D': 11,
        'H': 12,
        'V': 13,
        'N': 14,
        'X': 15
    }

struct2int = {
    'f': 0,
    't': 1,
    'i': 2,
    'h': 3,
    'm': 4,
    's': 5,
    'X': 6
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

    if not os.path.exists(os.path.join(dirname, 'dataset', type)):
        os.makedirs(os.path.join(dirname, 'dataset', type))

    if path_to_structures is None:
        path_to_structures = os.path.join(dirname, 'dataset', type, 'secondary_structure.json')

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
            row['structure'] = row['structure'].lower().ljust(max_seq_len, 'X') 

            # One hot encoding of the sequence, with padding as 'X' 
            for j, base in enumerate(row['sequence']):
                sequences[i,j] = seq2int[base]

            # One hot encoding of the structure, with padding as 'X' 
            for j, struct in enumerate(row['structure']):
                structures[i,j] = struct2int[struct]

        # Save sequences and structures as numpy arrays
        if save:
            np.save(save_path[0], sequences)
            np.save(save_path[1], structures)

        
        return sequences, structures



if __name__ == '__main__':

    sequences, structures = import_structure(type='full', save=True, reload=False)
    print("Loaded full dataset with shape: \n", sequences.shape, "\n", structures.shape)

    # sequences, structures = import_structure(type='train', save=True, reload=False)
    # print("Loaded train dataset with shape: \n", sequences.shape, "\n", structures.shape)

    # sequences, structures = import_structure(type='test', save=True, reload=False)
    # print("Loaded test dataset with shape: \n", sequences.shape, "\n", structures.shape)
