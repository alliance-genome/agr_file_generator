import sys, os
import wget
#import shutil
from pyfaidx import Fasta

class AssemblySequence(object):
    def __init__(self, assembly):
        self.assembly = assembly

        local_assembly_dir = "sequences/"
        local_filepath = local_assembly_dir +  assembly + ".fa"

        if not os.path.isfile(local_filepath):
            print(local_filepath)
            print("file does not exists")
            AssemblySequence.__downloadAssembly(assembly, local_assembly_dir, local_filepath)

        self.fa = Fasta(local_filepath)

    def __downloadAssembly(assembly, local_assembly_dir, local_filepath):
        assembly_to_s3_dict = {
                "R6.27": "https://s3.amazonaws.com/agrjbrowse/FlyBase/fruitfly/fasta/",
                "GRCz11": "https://s3.amazonaws.com/agrjbrowse/zfin/zebrafish/fasta/",
                "WBcel235": "https://s3.amazonaws.com/agrjbrowse/WormBase/c_elegans_PRJNA13758/fasta/",
                "Rnor_6.0": "https://s3.amazonaws.com/agrjbrowse/RGD/rat/fasta/",
                "GRCm38": "https://s3.amazonaws.com/agrjbrowse/MGI/mouse/fasta/"}
        url = assembly_to_s3_dict[assembly] + assembly + ".fa"
        if not os.path.exists(local_assembly_dir):
            os.makedirs(local_assembly_dir)
        wget.download(url, out=local_assembly_dir)

    def get(self, chromosome, start, end):
        start = start - 1 
        if chromosome.startswith("chr"):
            chromosome_str = chromosome[3:]
        else:
            chromosome_str = chromosome
        return self.fa[chromosome_str][start:end].seq
