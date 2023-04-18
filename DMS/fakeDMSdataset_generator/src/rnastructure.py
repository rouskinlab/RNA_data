import os
import numpy as np
import pandas as pd

def run_command(cmd):
    import subprocess
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return output.decode('utf-8')

class RNAstructure(object): 
    def __init__(self, rnastructure_path) -> None:
        self.rnastructure_path = os.path.join(rnastructure_path)

        dir_name = os.path.dirname(os.path.abspath(__file__))
        self.directory = os.path.join(dir_name, '..', 'temp')

    def predict_partition(self, temperature_k =None):
        # predict the partition of rna structures
        cmd = f"{os.path.join(self.rnastructure_path, 'partition')} {self.fasta_file} {self.pfs_file}"
        if temperature_k != None:
            cmd += ' --temperature '+str(temperature_k)
        run_command(cmd)
        
        # sum it into pairing probability
        run_command(os.path.join(self.rnastructure_path,'ProbabilityPlot')+ ' ' +self.pfs_file + ' -t '+self.prob_file)
        
        # read the probability file 
        with open(self.prob_file,"r") as f:
            lines=f.readlines()
            pairingPrediction={'i':[],'j':[],'p':[]}
            for x in range(len(lines)):
                if x>1:
                    ls=lines[x].split("\t")
                    pairingPrediction["i"]+=[int(ls[0])]
                    pairingPrediction["j"]+=[int(ls[1])]
                    pairingPrediction["p"]+=[float(ls[2])]
        df = pd.DataFrame(pairingPrediction)

        # turn log probability into probability
        df['p'] = df['p'].apply(lambda x: np.power(10, -x))
        
        # add the reverse complement so that when we sum, we get the probability of each base being paired with anyone
        df = pd.concat([df, df.rename(columns={'i':'j', 'j':'i'})])
        
        # cast: two bases being paired together => probability for each base of being paired with anyone 
        df = df.groupby(['i']).sum().drop(columns=['j'])

        # add p=0 for bases that are not paired with anyone
        df = df.reindex(range(1, len(self.sequence)+1), fill_value=0)

        # Logs sum can lead to a probability > 1, but it should never be > 1.05. So we'll sanity check that, then cap it at 1.
        assert (df['p'] >= 0).all() and (df['p'] <= 1.05).all(), 'The probability is not between 0 and 1.05, something is wrong. Check the log sum.'
        df['p'] = df['p'].apply(lambda x: min(x, 1))
        
        # Make sure all the bases are here
        assert len(df) == len(self.sequence), 'The number of bases in the sequence is not the same as the number of bases in the prediction.'
        
        return df['p'].tolist()
        
    def __make_temp_folder(self):
        isExist = os.path.exists(self.directory)
        if not isExist:
            os.makedirs(self.directory)
        return self.directory

    def __make_files(self, temp_prefix='temp'):
        self.pfs_file = f"{self.directory}/{temp_prefix}.pfs"
        self.ct_file = f"{self.directory}/{temp_prefix}.ct"
        self.dot_file = f"{self.directory}/{temp_prefix}_dot.txt"
        self.fasta_file = self.directory+'/'+temp_prefix+'.fasta'
        self.prob_file = self.directory+'/'+temp_prefix+'_prob.txt'

    def __create_fasta_file(self, reference, sequence):
        # push the ref into a temp file
        temp_fasta = open(self.fasta_file, 'w')
        temp_fasta.write('>'+reference+'\n'+sequence)
        temp_fasta.close()

    def predictPairingProbability(self, sequence, reference = 'reference'):
        self.sequence = sequence
        self.__make_temp_folder()
        self.__make_files()
        self.__create_fasta_file(reference, sequence)
        return self.predict_partition()
    

if __name__ == "__main__":
    rna = RNAstructure('/Users/ymdt/src/RNAstructure/exe')
    print("One line command:", rna.predictPairingProbability('AAGATATTCGAAAGAATATCTT'))

    