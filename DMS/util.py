import pandas as pd
import os
import numpy as np
import json
import sys

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

def import_DMS(path_to_dms, size=None, data_type='fake-data', save=False, reload=True):

    assert data_type in ['fake-data', 'sarah-testset'], "Type must be 'fake-data' or 'sarah-testset'"

    # Paths to the dataset
    dirname = os.path.dirname(os.path.abspath(__file__))
    save_path = [os.path.join(dirname, 'dataset', data_type, 'processed_sequences.npy'),
                 os.path.join(dirname, 'dataset', data_type, 'processed_dms.npy')]

    if path_to_dms is None:
        path_to_dms = os.path.join(dirname, 'dataset', data_type, 'dms_signal.json')

    # Import the dataset and check size
    df = pd.read_json(path_to_dms).T
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
            dms_signals = np.load(save_path[1])

            if len(sequences) < size or len(dms_signals) < size:
                print("Dataset too small, creating new one")
                return import_DMS(path_to_dms, data_type=data_type, size=size, save=save, reload=False)
            else:
                idx = np.random.choice(len(sequences), size=size, replace=False)
                if save:
                    np.save(save_path[0], sequences[idx])
                    np.save(save_path[1], dms_signals[idx])
                
                return sequences[idx], dms_signals[idx]
        else:
            print("Dataset not found, creating new one")
            return import_DMS(path_to_dms, data_type=data_type, size=size, save=save, reload=False)
        

    else:        
        
        # Get a random sample of the dataset
        idx = np.random.choice(len(df), size=size, replace=False)
        df = df.iloc[idx].reset_index(drop=True)

        # Get max sequence length
        max_seq_len = df['sequence'].str.len().max()

        # Init numpy arrays for sequences and structures
        sequences = np.zeros((len(df), max_seq_len))
        dms_signals = np.zeros((len(df), max_seq_len))

        # Iterate over the rows of the dataframe
        for i, row in df.iterrows():

            # Print progress
            if i%1000:
                sys.stdout.write("Processing dataset: %d%%   \r" % (100*i/len(df)) )
                sys.stdout.flush()
            
            # Apply zero padding to each row (sequence and dms_signal)
            row['sequence'] = row['sequence'].upper().ljust(max_seq_len, 'X')
            dms_signals[i, :] = row['dms_signal'] + [0]*(max_seq_len-len(row['dms_signal']))
            
            # Integer encoding of the sequence, with padding as 'X' 
            for j, base in enumerate(row['sequence']):
                sequences[i,j] = seq2int[base]

        # Save sequences and structures as numpy arrays
        if save:
            np.save(save_path[0], sequences)
            np.save(save_path[1], dms_signals)

        
        return sequences, dms_signals


if __name__ == '__main__':

    sequences, dms_signals = import_DMS('/Users/alberic/Desktop/Pro/RouskinLab/projects/deep_learning/RNA_data/DMS/dataset/fake-data/dms_signal.json', 
                    data_type='fake-data', size=100000, save=True, reload=False)
    
    print("Loaded full dataset with shape: \n", sequences.shape, "\n", dms_signals.shape)
