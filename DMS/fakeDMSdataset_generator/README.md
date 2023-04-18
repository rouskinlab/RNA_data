# fakeDMSdata

Generate fake DMS signal using [RNAstructure pairing probabilities](https://rna.urmc.rochester.edu/Text/ProbabilityPlot.html), to train our models until we get the experimental data.

The function `createFakeData` reads a fasta file and outputs a dict with reference, sequence and synthetic DMS signal. It is located in `createFakeData.py`.

## Function description

```Python
def createFakeData(fasta_file, rnastructure_path='', sequencer_noise=0.001):
    """
    Reads a fasta file and outputs a dict with reference, sequence and synthetic DMS signal.

    Args:
        fasta_file (str): The path to the input fasta file.
        rnastructure_path (str): The path to the RNAstructure executable. Default is '' (i.e. RNAstructure is in the PATH).
        sequencer_noise (float): The amount of sequencer noise to add to the DMS signal. Default is 0.001.
    
    Returns:
        data (dict): A dictionary with reference as key and sequence / synthetic DMS signal as value.
    
    Example:
        >>> createFakeData('testData/refs.fasta')
        Predicting RNA structures: 100%|███████████████████████████████████| 3/3 [00:01<00:00,  1.92seq/s]
        {'3042-O-flank_1=hp1-DB': 
            {
            'sequence': 'TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCAAGATATTCGAAAGAATATCTTTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTGTATTACGAGTTCGCTACTCGTTCCTTTCGA', 
            'fakeDMS': [0.0008763192815798796, 0.003282160575810705, 0.019867824327687585, 0.02829582753478791, 0.031266714341378524, 0.010288700618344545, 0.22755049833180016, 0.5741991237121613, 0.5018175337184233, 0.004588974386334777, 0.0027875607249180545, 0.010719083819803172, 0.0311800339547883, 0.026734132953685468, 0.5171140931132598, 0.526320420155328, 0.5302629606137425, 0.5973646221995781, 0.06787639095534494, 0.996380706717648, 0.9988325172245213, 0.9986708144405785, 0.9984979885954683, 0.9984931638770012, 0.9984703346038281, 0.9986003515173016, 0.0005489978550856028,
            },
        '3043-CC-flank_1=hp1-DB_2':
            ...
        
    """

```
