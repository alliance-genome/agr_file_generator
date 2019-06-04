import os

import wget
from pyfaidx import Fasta

class AssemblySequence:

    def __init__(self, assembly, local_assembly_dir='sequences'):
        self.assembly = assembly
        local_filepath = os.path.join(local_assembly_dir, assembly + '.fa')
        if not os.path.isfile(local_filepath):
            print('Assembly file does not exist:')
            self.download(assembly, local_assembly_dir, local_filepath)
        self.fa = Fasta(local_filepath)

    @classmethod
    def download(cls, assembly, local_assembly_dir, local_filepath):
        assembly_to_s3_dict = {
            'R6.27': 'https://s3.amazonaws.com/agrjbrowse/FlyBase/fruitfly/fasta/',
            'GRCz11': 'https://s3.amazonaws.com/agrjbrowse/zfin/zebrafish/fasta/',
            'WBcel235': 'https://s3.amazonaws.com/agrjbrowse/WormBase/c_elegans_PRJNA13758/fasta/',
            'Rnor_6.0': 'https://s3.amazonaws.com/agrjbrowse/RGD/rat/fasta/',
            'GRCm38': 'https://s3.amazonaws.com/agrjbrowse/MGI/mouse/fasta/'
        }
        url = assembly_to_s3_dict[assembly] + assembly + '.fa'
        if not os.path.exists(local_assembly_dir):
            os.makedirs(local_assembly_dir)
        wget.download(url, out=local_assembly_dir)

    def get(self, chromosome, start, end):
        start -= 1
        if chromosome.startswith('chr'):
            chromosome_str = chromosome[3:]
        else:
            chromosome_str = chromosome
        return self.fa[chromosome_str][start:end].seq
