import os, json
from .rnastructure import RNAstructure

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

        return Ct.get_reference_from_title(ct_file), sequence.upper().replace('T', 'U'), paired_bases

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
    def parse(dreem_output_file):
        """Parse a dreem output file and return the references and sequences"""
        data = json.load(open(dreem_output_file, 'r'))
        for ref, v in data.items():
            if type(v) != dict:
                continue
            yield ref, v['full']['sequence'], v['full']['pop_avg']['sub_rate']

