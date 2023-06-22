import os
import numpy as np
import sys
import pandas as pd
import json
from tqdm import tqdm

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'RNAstructure'))
from src.util import fastaToDict
from util import seq2int
from RNAstructure import predictFromSequence



def predict_pairing_matrix(input_dir=None, path_to_RNAStructure='', constraint=False):

    # Read the sample with DMS constraints
    if constraint:
        with open(input_dir) as json_file:
            data = json.load(json_file)
    
    # Or only use the sequences from fasta file
    else: 
        data = fastaToDict(input_dir)
    

    # Go over each sequence and predict its structure
    ref_todelete = []
    for ref in tqdm(data.keys()):

        # Check if the sequence contains only A, C, G, U/T, discard it from the dataset otherwise
        bases = ['A','U','C','G','T','N','a','u','c','g','t','n']
        if not all(c in bases for c in data[ref]['sequence']):
            print(f'WARNING: sequence {ref} contains non-canonical bases, skipping it')
            ref_todelete.append(ref)
            continue
        
        pairing_constraint = []
        if constraint:
            pairing_constraint = data[ref]['dms_signal']

        prediction = predictFromSequence(data[ref]['sequence'], dms = pairing_constraint,
                                        predict_pairs=True, predict_structure=False, predict_pairing_probability=False, 
                                        rnastructure_path=path_to_RNAStructure)


        paired_bases = np.array(prediction['list_pairs'])-1
        paired_bases = np.unique(np.sort(paired_bases, axis=1), axis=0)
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
        
    return pd.DataFrame.from_dict(data, orient='index')

def convert_CT_to_pairing_matrix(path_to_ct=None):

    if path_to_ct is None:
        path_to_ct = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dataset', 'test_PDB', 'CT_files_pdbee')

    # Read all CT files to get the list of sequences and pairing matrices
    data = {}
    for file in os.listdir(path_to_ct):
        if file.endswith(".txt") or file.endswith(".ct"):
            with open(os.path.join(path_to_ct, file)) as f:
                data_file = f.readlines()

                sequence = []
                base_pair = []
                for row in data_file[1:]:

                    if row=='':
                            continue
                    
                    row_data = row.split()
                    
                    if row_data[1] == 'T' or row_data[1] == 't':
                        sequence.append('U') # replace T by U
                    else:
                        sequence.append(row_data[1].upper())
                    if int(row_data[4]) != 0:
                        base_pair.append(sorted([int(row_data[0])-1, int(row_data[4])-1])) # Sort and -1 to get 0-based indexing

                # Check if the sequneces and structures are valid
                base_pair = np.unique(base_pair, axis=0)
                if len(base_pair) > 0:
                    if (all((b in seq2int.keys()) and (b!='X') for b in sequence) 
                        and (base_pair<len(sequence)).all() 
                        and len(np.unique(base_pair[:,0]))==len(base_pair)
                        and len(np.unique(base_pair[:,1]))==len(base_pair) ):
                        data[file.split('.')[0]] = {'sequence': ''.join(sequence), 'paired_bases': base_pair.tolist()}
                else:
                    if (all((b in seq2int.keys()) and (b!='X') for b in sequence)):
                        data[file.split('.')[0]] = {'sequence': ''.join(sequence), 'paired_bases': []}


    # Postprocess data to remove duplicate sequences with different structures
    df = pd.DataFrame.from_dict(data, orient='index')

    unique_seqs = {}
    for i, seq in enumerate(df.sequence):

        struct = np.array(df['paired_bases'][i])

        if seq not in unique_seqs.keys():
            unique_seqs[seq] = [struct]
        else:
            unique_seqs[seq].append(struct)

    # Only keep sequence duplicate that have the same pairing matrix
    for seq in unique_seqs.keys():
        if len(unique_seqs[seq]) > 1:
            matches = []
            for i in range(len(unique_seqs[seq])):
                for j in range(i+1, len(unique_seqs[seq])):
                    matches.append(np.array_equal(unique_seqs[seq][i], unique_seqs[seq][j]))

            # Remove sequences from dataframe if f1 score is not 1.0, keep only one of the structures otherwise
            if not np.array(matches).all():
                df = df[df['sequence'] != seq]
            else:
                row = df[df['sequence'] == seq].head(1)
                df = df[df['sequence'] != seq]
                df = pd.concat([df, row])

    return df



# Generate secondary structure dataset
# Arguments:
#   - args[1]: dataset_type: predict or test_PDB
#   - args[2]: input_dir: path to the input directory. If None use the data inside this repository
#   - args[3]: output_dir: path to the output directory. If None save the data in dataset/input_dir
#   - args[4]: path_to_RNAStructure: path to the RNAStructure executable

if __name__ == '__main__':

    # Get args
    dataset_type = 'predict' 
    input_dir = '/Users/alberic/Desktop/Pro/RouskinLab/projects/deep_learning/RNA_data/DMS/dataset/sarah_supermodels/dms_signal.json'
    output_dir = '/Users/alberic/Desktop/Pro/RouskinLab/projects/deep_learning/RNA_data/secondary_structure/dataset/test_Sara'
    path_to_RNAStructure = '/Users/alberic/RNAstructure/exe/'

    # Get the secondary structure sequences as pairing matrix
    if dataset_type == 'predict':
        df = predict_pairing_matrix(input_dir=input_dir, path_to_RNAStructure=path_to_RNAStructure, constraint=True)
    elif dataset_type == 'test_PDB':
        df = convert_CT_to_pairing_matrix(path_to_ct=input_dir)
    else:
        raise ValueError('ERROR: dataset type not recognized')

    # Save the dataframes as json
    if output_dir is None:
        dir_name = os.path.dirname(os.path.abspath(__file__))
        output_dir = os.path.join(dir_name, 'dataset', dataset_type)
    
    # Write dataframe data to json file
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    df.to_json(os.path.join(output_dir, 'secondary_structure.json'), orient='index', indent=2)