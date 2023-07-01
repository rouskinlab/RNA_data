import os

# Define the one-hot encodings for the sequences and structures
seq2int = {
        'A': 1,
        'C': 2,
        'G': 3,
        'T': 4,
        'U': 4,
        'N': 5,
        'Y': 6,
        'R': 7,
        'K': 8,
        'W': 9,
        'S': 10,
        'M': 11,
        'B': 12,
        'D': 13,
        'H': 14,
        'V': 15,
        'X': 0
    }

dot2int = {'.': 1, '(': 2, ')': 3, 'X': 0}
int2dot = ['X', '.', '(', ')']


def fastaToDict(fasta_file):
    """
    Reads a fasta file and returns a dictionary {sequence: reference}.

    Args:
        fasta_file (str): The path to the input fasta file.
    Returns:
        refSeq (dict): A dictionary with reference as key and sequence as value.


    Example:
    >>> import os
    >>> os.makedirs('temp', exist_ok=True)
    >>> f = open('temp/test.fasta', 'w')
    >>> f.write('>ref1\\nACGT\\n>ref2\\nCGTA')
    21
    >>> f.close()
    >>> print(fastaToDict('temp/test.fasta'))
    {'ref1': {'sequence': 'ACGT'}, 'ref2': {'sequence': 'CGTA'}}

    """

    refSeq = {}

    assert os.path.exists(fasta_file), f'The fasta file {fasta_file} does not exist.'
    with open(fasta_file, 'r') as handle:
        while True: # while not EOF
            ref, seq = handle.readline(), handle.readline()
            if ref == '' or seq == '':
                break
            assert ref[0] == '>', 'The reference should start with ">"'

            refSeq[ref[1:].strip()] = {'sequence': seq.strip().upper()}

    return refSeq


def add_braces_if_no_braces(s:str):
    if not s[0] == '{' or not s[-1] == '}':
        return '{' + s + '}'
    return s


def source_env(path):
    """
    Source the environment variables from the file at path.

    Args:
        path (str): The path to the file to source.
    """
    out = {}
    with open(path, 'r') as f:
        for line in f.readlines():
            line = line.split('#')[0].strip()
            line = line.replace('export', '')
            key, value = line.split('=')
            key, value = key.replace(' ','').replace('"','').strip(), value.replace(' ','').replace('"','').strip()
            os.environ[key] = value
            out[key] = value

    return out