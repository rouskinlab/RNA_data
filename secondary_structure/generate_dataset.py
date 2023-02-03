import os
import numpy as np
import subprocess
import re
import sys
import lib_forgi
import pandas as pd

# Function from https://github.com/kkyamada/bert-rbp
def get_mea_structures(path_to_fasta, path_to_linearpartition, training=False):
    
    command = 'cat {} | {} -M'.format(path_to_fasta, path_to_linearpartition)
    
    # could use capture_output=True to not print to stdout, but we will need to print progress bar another way
    output = subprocess.run(command, shell=True, stdout=subprocess.PIPE, check=True)
    
    assert output.returncode == 0, 'ERROR: in execute_linearpartition_mea for {}'.format(path_to_fasta)

    # Convert outputs of linearpartition to integer encoding
    def make_node_set(numbers):
        numbers = list(map(int, numbers))
        ans = set()
        while len(numbers) > 1:
            a, b = numbers[:2]
            numbers = numbers[2:]
            for n in range(a - 1, b):
                ans.add(n)  # should be range a,b+1 but the biologists are weired
        return ans
    
    structures = output.stdout.decode().strip().split('\n')
    
    # f: 'dangling start', 't': 'dangling end', 'i': 'internal loop', 'h': 'hairpin loop', 'm': 'multi loop', 's': 'stem'
    entity_lookup = {'f': 0, 't': 1, 'i': 2, 'h': 3, 'm': 4, 's': 5}

    structure_seqs = []
    sequences = []
    for i_s, structure in enumerate(structures):

        # If the structure is the sequence information:
        if i_s%4 == 0:
            sequences.append(structure)

        # If the structure is a dot-bracket structure:
        if (i_s-2)%4 == 0:
            bg = lib_forgi.BulgeGraph()
            bg.from_dotbracket(structure, None)
            forgi = bg.to_bg_string()

            # Save the structure sequence as a numpy array of characters
            structure_sequence = [0]*len(structure)
            for line in forgi.split('\n')[:-1]:
                # if the line starts with 'define' we know that annotation follows...
                if line[0] == 'd':
                    l = line.split()
                    # first we see the type
                    entity = l[1][0]
                    # then we see a list of nodes of that type.
                    entity_index = entity_lookup[entity]
                    for n in make_node_set(l[2:]):
                        structure_sequence[n] = entity

            structure_sequence = ''.join(structure_sequence)
            structure_seqs.append(structure_sequence)
                        
    return sequences, structure_seqs

if __name__ == '__main__':

    dir_name = os.path.dirname(os.path.abspath(__file__))

    # Get path to fasta file
    if len(sys.argv) >= 2:
        path_to_fasta = sys.argv[1]
    else:
        path_to_fasta = os.path.join(dir_name, '..', 'dataset', 'sequences_small.fasta')

    # Compile LinearPartition if it hasn't been compiled yet
    if not os.path.exists(os.path.join(dir_name, 'LinearPartition', 'bin')):
        print('Compiling LinearPartition')
        os.chdir(os.path.join(dir_name, 'LinearPartition'))
        os.system('make')
        os.chdir(dir_name)

    # Get the secondary structure sequences
    path_to_linearpartition = os.path.join(dir_name, 'LinearPartition', 'linearpartition')
    sequences, structure_seqs = get_mea_structures(path_to_fasta, path_to_linearpartition, training=False)
    
    # Save the secondary structure sequences and sequences as json
    df = pd.DataFrame({'sequence': sequences, 'structure': structure_seqs})
    df.to_json(os.path.join(dir_name, 'dataset', 'secondary_structure.json'), indent=2)