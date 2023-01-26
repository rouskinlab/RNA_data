import os
import numpy as np
import subprocess
import re
import sys
import lib_forgi

# Function from https://github.com/kkyamada/bert-rbp
def get_mea_structures(path_to_fasta, path_to_linearpartition, training=False):
    # path_to_fasta = os.path.join(path_to_file, "dev_bpp.fasta")

    
    if not 'linearpartition' in path_to_linearpartition:
        path_to_linearpartition = os.path.join(path_to_linearpartition, "linearpartition")
    
    #if os.path.isfile(path_to_mea):
    #    command = 'rm {}'.format(path_to_mea)
    #    subprocess.run(command, shell=True)
    # command = 'cat {} | {} -m > {}'.format(path_to_fasta, path_to_linearpartition, path_to_mea)
    command = 'cat {} | {} -M'.format(path_to_fasta, path_to_linearpartition)
    # print(command)
    #process_result = 1
    #with open(path_to_mea, 'w') as output_file:
    #    process_result = subprocess.run(command, shell=True, stdout=output_file)
    # print(command)
    output = subprocess.run(command, shell=True, stdout=subprocess.PIPE, check=True)
    # process_result = subprocess.run(command, shell=True)
    # print("process done")
    
    # structure_seqs = np.array([])
    structure_seqs = []
    if output.returncode == 1:
        print('ERROR: in execute_linearpartition_mea for {}'.format(path_to_fasta))
    else:
        def make_node_set(numbers):
            numbers = list(map(int, numbers))
            ans = set()
            while len(numbers) > 1:
                a, b = numbers[:2]
                numbers = numbers[2:]
                for n in range(a - 1, b):
                    ans.add(n)  # should be range a,b+1 but the biologists are weired
            return ans
        
        #structures = []
        #with open(path_to_mea, 'r') as f:
        #    structures = f.readlines()
        structures = output.stdout.decode().strip().split('\n')
        
        count = 0
        entity_lookup = {'f': 0, 't': 1, 'i': 2, 'h': 3, 'm': 4, 's': 5}
        # f: 'dangling start', 't': 'dangling end', 'i': 'internal loop', 'h': 'hairpin loop', 'm': 'multi loop', 's': 'stem'

        for structure in structures:
            #structure = structure.strip()
            if re.fullmatch("[\.\(\)]+", structure):
                #print(structure)
                bg = lib_forgi.BulgeGraph()
                bg.from_dotbracket(structure, None)
                forgi = bg.to_bg_string()
                structure_sequence = np.zeros([6, len(structure)])
                for line in forgi.split('\n')[:-1]:
                    # if the line starts with 'define' we know that annotation follows...
                    if line[0] == 'd':
                        l = line.split()
                        # first we see the type
                        entity = l[1][0]
                        # then we see a list of nodes of that type.
                        entity_index = entity_lookup[entity]
                        for n in make_node_set(l[2:]):
                            structure_sequence[entity_index,n] = 1
                structure_seqs.append(structure_sequence)
                            
                # if len(structure_seqs)==0:
                #     structure_seqs = structure_sequence
                # else:
                #     structure_seqs = np.concatenate([structure_seqs, structure_sequence])
    return structure_seqs

if __name__ == '__main__':
    # path_to_fasta = os.path.join(path_to_file, "dev_bpp.fasta")
    path_to_fasta = sys.argv[1]
    path_to_linearpartition = sys.argv[2]
    structure_seqs = get_mea_structures(path_to_fasta, path_to_linearpartition, training=False)
    
    with open("datasets/secondary_structure.txt", "ab") as f:
        for seq in structure_seqs:
            f.write(b"\n")
            np.savetxt(f, seq)