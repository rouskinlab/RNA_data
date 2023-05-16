import pandas as pd
import os
import numpy as np

import sys
import torch

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

dot2int = {'.': 1, '(': 2, ')': 3, 'X': 0}
int2dot = ['X', '.', '(', ')']

def import_structure(path_to_structures=None, size=None, save=False, reload=True):
    """
    Import the secondary structure dataset and convert to pairing matrix, with padding.
    Each pairing matrix is directly saved as a separate npy file.

    Each row of the dataset contains a sequence and a structure.
    The sequences contains the nucleotides A, C, G, T, U, N
    The structures is binary pairing matrix with 0 for unpaired bases and 1 for paired bases

    :param path_to_structures: Path to the secondary structure dataset in json format
    :param size: The number of datapoints to import
    :param save: Whether to save the dataset as a numpy array
    :param reload: Whether to reload the dataset from the numpy array

    :return: The sequences as numpy array and the list of numpy file names of the pairing matrices
    """

    # Paths to the dataset
    dirname = os.path.dirname(os.path.abspath(__file__))
    save_path = [os.path.join(dirname, 'dataset', 'pairing_matrix', 'processed_sequences.npy'),
                 os.path.join(dirname, 'dataset', 'pairing_matrix', 'processed_structures'),
                 os.path.join(dirname, 'dataset', 'pairing_matrix', 'processed_structures', 'references.txt')]
    
    if not os.path.exists(save_path[1]):
        os.makedirs(save_path[1])

    if path_to_structures is None:
        path_to_structures = os.path.join(dirname, 'dataset', 'pairing_matrix', 'secondary_structure.json')


    # Check if the dataset is already saved
    if reload:
        
        # Import the dataset and check size
        if size is None:
            df = pd.read_json(path_to_structures).T
            size = len(df)

        if os.path.exists(save_path[0]) and os.path.exists(save_path[2]):
            print("Loading saved dataset")
            sequences = np.load(save_path[0])
            references = np.loadtxt(save_path[2], dtype=str)

            if len(sequences) < size or len(references) < size:
                print("Dataset too small, creating new one")
                return import_structure(path_to_structures, size=size, save=save, reload=False)
            else:
                idx = np.random.choice(len(sequences), size=size, replace=False)
                if save:
                    np.save(save_path[0], sequences[idx])
                    np.save(save_path[2], references[idx])
                
                return sequences[idx], references[idx]
        else:
            print("Dataset not found, creating new one")
            return import_structure(path_to_structures, size=size, save=save, reload=False)
        

    else:  

        # Import the dataset and check size
        df = pd.read_json(path_to_structures).T
        if size is None:
            size = len(df)      
        elif size > len(df):
            print("Requested size too large, using full dataset")
            size = len(df)
        
        # Get a random sample of the dataset
        idx = np.random.choice(len(df), size=size, replace=False)
        df = df.iloc[idx]

        # Get max sequence length
        max_seq_len = df['sequence'].str.len().max()

        # Init numpy arrays for sequences
        sequences = np.zeros((len(df), max_seq_len))
        references = []

        # Iterate over the rows of the dataframe
        for i, (reference, row) in enumerate(df.iterrows()):

            # Print progress
            if i%1000:
                sys.stdout.write("Processing dataset: %d%%   \r" % (100*i/len(df)) )
                sys.stdout.flush()
            
            # Apply zero padding to each row (sequence and structure)
            row['sequence'] = row['sequence'].upper().ljust(max_seq_len, 'X')

            # Integer encoding of the sequence, with padding as 'X' 
            for j, base in enumerate(row['sequence']):
                sequences[i,j] = seq2int[base]

            # Create pairing matrix from list of paired bases
            paired_bases = np.array(row['paired_bases'])
            pairing_matrix = torch.zeros((max_seq_len, max_seq_len)).bool()
            if len(paired_bases) > 0:
                pairing_matrix[paired_bases[:,0], paired_bases[:,1]] = 1.0
                pairing_matrix[paired_bases[:,1], paired_bases[:,0]] = 1.0

            references.append(reference)
            if save:
                torch.save(pairing_matrix, os.path.join(save_path[1], reference+'.pt'))

        # Save sequences and structures as numpy arrays
        if save:
            np.save(save_path[0], sequences)
            np.savetxt(save_path[2], references, fmt='%s')

        
        return sequences, references



if __name__ == '__main__':

    sequences, references = import_structure(save=True, reload=False)
    print("Loaded full dataset with shape: \n", sequences.shape, "\n", len(references))