import numpy as np
import os 

# Read the fasta file
dir_name = os.path.dirname(os.path.abspath(__file__))
path_to_fasta = os.path.join(dir_name, 'sequences_half.fasta')
sequences = []
names = []
with open(path_to_fasta, 'r') as f:
    lines = f.readlines()
    for line in lines:
        if line[0] != '>':
            sequences.append(line.strip().upper())
        else:
            names.append(line.strip())


print('Number of sequences: {}'.format(len(sequences)))
print('Number of unique sequences: {}'.format(len(np.unique(sequences))))

print('Removing duplicate sequences')
_, idx = np.unique(sequences, return_index=True)
idx = sorted(idx)

sequences = np.array(sequences)[idx]
names = np.array(names)[idx]

print("Removing duplicate names")
_, idx = np.unique(names, return_index=True)
idx = sorted(idx)

sequences = sequences[idx]
names = names[idx]

print('Number of sequences: {}'.format(len(sequences)))
print('Number of unique sequences: {}'.format(len(np.unique(sequences))))

print('Saving sequences')
with open(path_to_fasta, 'w') as f:
    for name, seq in zip(names, sequences):
        f.write(name+'\n')
        f.write(seq+'\n')