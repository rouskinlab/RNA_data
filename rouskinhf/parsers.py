import os, json
import numpy as np
from .rnastructure import RNAstructure
from .util import DreemUtils
import pandas as pd


class Ct:
    def parse(ct_file):
        """Parse a ct file and return the sequence and structure

        Args:
            ct_file (str): path to ct file

        Returns:
            (str,str,str): (reference, sequence, paired_bases)


        """
        with open(ct_file, 'r') as f:
            lines = f.readlines()

        paired_bases, sequence = [], ''
        for line in lines[1:]:
            if line.strip() == '':
                break
            utr5, base, _, _, utr3, _ = line.split()
            sequence += base
            if int(utr3) > int(utr5) and int(utr3) != 0:
                paired_bases.append([int(utr5)-1, int(utr3)-1])

        return Ct.get_reference_from_title(ct_file), sequence, paired_bases


    def parse_list(ct_files):
        """Parse a list of ct files and return the sequences and structures"""
        return [Ct.parse(ct_file) for ct_file in ct_files]

    def get_reference_from_title(ct_file):
        return os.path.basename(ct_file).split('.')[0]

    def predict_dms(ct_file):
        """Predict the dms of a ct file using RNAstructure"""
        rna = RNAstructure()
        return rna.predictPairingProbability(Ct.parse(ct_file)[1])


class Fasta:
    def parse(fasta_file):
        """Parse a fasta file and return the references and sequences"""
        refs, seqs = [], []
        with open(fasta_file, 'r') as f:
            for line in f.readlines():
                if line[0] == '>':
                    refs.append(line[1:].strip())
                else:
                    seqs.append(line.strip())
        assert len(refs) == len(seqs), "The number of references and sequences in the fasta file must be the same ({} vs {})".format(len(refs), len(seqs))
        return seqs, refs

    def get_name(fasta_file):
        return os.path.basename(fasta_file).split('.')[0]

    def predict_structure(sequence):
        """Predict the structure of a sequence using RNAstructure"""
        rna = RNAstructure()
        return rna.predictStructure(sequence)

    def predict_dms(sequence):
        """Predict the dms of a sequence using RNAstructure"""
        rna = RNAstructure()
        return rna.predictPairingProbability(sequence)


class DreemOutput:
    def parse(dreem_output_file, drop_duplicates=True, max_mut=0.5):
        """Parse a dreem output file and return the references and sequences.

        Args:
            dreem_output_file (str): path to dreem output file
            drop_duplicates (bool): drop duplicate sequences
            max_mut (float): drop sequences with at least one mutation > max_mut

        Returns:
            (str,str,str): (reference, sequence, sub_rate)
        """

        df = pd.DataFrame(DreemUtils.flatten_json(DreemUtils.sort_dict(json.load(open(dreem_output_file, 'r')))))[['reference', 'sequence', 'sub_rate']]

        df['max_mut'] = df.sub_rate.apply(lambda x: max(x))
        len_before = len(df)
        df = df[df.max_mut < max_mut].reset_index(drop=True)

        len_before = len(df)
        if drop_duplicates:
            df.drop_duplicates(subset='sequence', inplace=True)


        values = np.concatenate(df.sub_rate.values)
        percentile975 = np.percentile(values, 97.5)

        def normalize(sub_rate):
            sub_rate = np.array(sub_rate)
            sub_rate = sub_rate / percentile975
            sub_rate[sub_rate > 1] = 1
            return sub_rate.tolist()

        for _, row in df.iterrows():
            yield row['reference'], row['sequence'], normalize(row['sub_rate'])

