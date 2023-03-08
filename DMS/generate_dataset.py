import os
import json
import dreem
import numpy as np
import pandas as pd

from plotly.subplots import make_subplots
import plotly.graph_objects as go

seq2int = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3,
        'U': 3,
        'Y': 4,
        'R': 5,
        'K': 6,
        'W': 7,
        'S': 8,
        'M': 9,
        'B': 10,
        'D': 11,
        'H': 12,
        'V': 13,
        'N': 14,
        'X': 15
    }


def get_dms_signals(path_data):

    data = []
    for f in os.listdir(path_data):
        if f.endswith('.json'):
            print(' ' + f, end=' ')
            data.append(json.load(open(os.path.join(path_data, f), 'r')))

    # Create a study object to process the data
    study = dreem.draw.study.Study(
        data = data
    )
    study.df = study.df[study.df['section']=='full']
    study.df = study.df[study.df['DMS_conc_mM']>0.0]

    # Filter out low quality data
    min_base_coverage = 1000
    study.df = study.df[study.df['min_cov'] > min_base_coverage].reset_index(drop=True)

    sequences=[]
    signals=[]

    # Normalize each sample
    for i, sample in enumerate(study.df['sample'].unique()):

        sub_df = study.get_df(sample=sample)

        # Concatenate all the DMS to compute the cutoff
        new_dms = np.concatenate(sub_df['sub_rate'].values) 
        max_dms = np.median(new_dms[new_dms>np.percentile(new_dms, 95)])
        
        # Save the values
        for seq, dms in zip( sub_df['sequence'], sub_df['sub_rate']):

            sequences.append([seq2int[base] for base in seq ])

            dms_norm = np.array(dms)/max_dms
            dms_norm[dms_norm>1] = 1
            signals.append(dms_norm)

    return sequences, signals


if __name__ == '__main__':

    # Get json output from DREEM
    dir_name = os.path.dirname(os.path.abspath(__file__))
    path_data = os.path.join(dir_name, 'dataset', 'samples')

    sequences, signals = get_dms_signals(path_data)

    df = pd.DataFrame({'sequence': sequences, 'DMS': signals})
    df.to_json(os.path.join(dir_name, 'dataset', 'dms_signals.json'), indent=2)


