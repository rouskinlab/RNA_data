def fastaToDict(fasta_file):
    """
    Reads a fasta file and returns a dictionary {sequence: reference}.
    
    Args:
        fasta_file (str): The path to the input fasta file.
    Returns:
        refSeq (dict): A dictionary with reference as key and sequence as value.
        
    """

    refSeq = {}

    with open(fasta_file, 'r') as handle:
        ref = 'placeholder'
        while True: # while not EOF
            ref, seq = handle.readline(), handle.readline()
            if ref == '' or seq == '':
                break
            assert ref[0] == '>', 'The reference should start with ">"'
            
            refSeq[ref[1:].strip()] = {'sequence': seq.strip().upper()}
    
    return refSeq