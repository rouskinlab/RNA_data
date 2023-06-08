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

def parse_args():
    # Dataset type: synthetic or test_PDB
    if len(sys.argv) >= 2:
        dataset_type = sys.argv[1]
    else:
        dataset_type = 'test_PDB'

    # Input directory
    if len(sys.argv) >= 3:
        input_dir = sys.argv[1]
    else:
        input_dir = None

    # Output directory
    if len(sys.argv) >= 4:
        output_dir = sys.argv[2]
    else:
        output_dir = None
    
    # Path to RNAStructure
    if len(sys.argv) >= 5:
        path_to_RNAStructure = sys.argv[3]
    else:
        path_to_RNAStructure = ''

    return dataset_type, input_dir, output_dir, path_to_RNAStructure


def predict_pairing_matrix(path_to_fasta=None, path_to_RNAStructure=''):

    if path_to_fasta is None:
        path_to_fasta = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sequence_dataset', 'sequences_half.fasta')

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
        
    return pd.DataFrame.from_dict(data, orient='index')

def convert_CT_to_pairing_matrix(path_to_ct=None):

    if path_to_ct is None:
        path_to_ct = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dataset', 'test_PDB', 'CT_files_pdbee')

    # Read all CT files to get the list of sequences and pairing matrices
    data = {}
    for file in os.listdir(path_to_ct):
        if file.endswith(".txt"):
            with open(os.path.join(path_to_ct, file)) as f:
                data_file = f.readlines()

                sequence = []
                base_pair = []
                for row in data_file[1:]:
                    row_data = row.split(' ')
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
#   - args[1]: dataset_type: synthetic or test_PDB
#   - args[2]: input_dir: path to the input directory. If None use the data inside this repository
#   - args[3]: output_dir: path to the output directory. If None save the data in dataset/input_dir
#   - args[4]: path_to_RNAStructure: path to the RNAStructure executable

if __name__ == '__main__':

    # Get args
    dataset_type, input_dir, output_dir, path_to_RNAStructure  = parse_args()

    # Get the secondary structure sequences as pairing matrix
    if dataset_type == 'synthetic':
        df = predict_pairing_matrix(path_to_fasta=input_dir, path_to_RNAStructure='/Users/alberic/RNAstructure/exe/')
    elif dataset_type == 'test_PDB':
        df = convert_CT_to_pairing_matrix(path_to_ct=input_dir)
    else:
        raise ValueError('ERROR: dataset type not recognized')

    # Save the dataframes as json
    if output_dir is None:
        dir_name = os.path.dirname(os.path.abspath(__file__))
        output_dir = os.path.join(dir_name, 'dataset', dataset_type)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    # Write dataframe data to json file
    df.to_json(os.path.join(output_dir, 'secondary_structure.json'), orient='index', indent=2)