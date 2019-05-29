import sys, os
import wget
import gzip
import shutil
from pyfaidx import Fasta

class AssemblySequence(object):
    def __init__(self, assembly, chromosome):
        self.assembly = assembly
        self.chromosome = chromosome

        if assembly not in ["GRCm38"]:
            assembly_zipped = True
        else:
            assembly_zipped = False

        local_assembly_dir = "sequences/" + assembly + "/"
        local_filepath = local_assembly_dir + chromosome + ".fa"

        if not os.path.isfile(local_filepath):
            print(local_filepath)
            print("file dos not exists")
            AssemblySequence.__downloadAssembly(assembly, chromosome, assembly_zipped, local_assembly_dir, local_filepath)

        self.fa = Fasta(local_filepath)

    def __downloadAssembly(assembly, chromosome, assembly_zipped, local_assembly_dir, local_filepath):
        assembly_to_s3_dict = {
                "R6.27": "https://s3.amazonaws.com/agrjbrowse/FlyBase/fruitfly/fasta/",
                "GRCz11": "https://s3.amazonaws.com/agrjbrowse/zfin/zebrafish/fasta/",
                "WBcel235": "https://s3.amazonaws.com/agrjbrowse/WormBase/c_elegans_PRJNA13758/fasta/",
                "Rnor_6.0": "https://s3.amazonaws.com/agrjbrowse/RGD/rat/fasta/",
                "GRCm38": "https://s3.amazonaws.com/agrjbrowse/MGI/mouse/fasta/"}
        if chromosome.startswith("chr"):
            chromosome_str = chromosome[3:]
            removed_chr = True
        else:
            chromosome_str = chromosome
            removed_chr = False
        url = assembly_to_s3_dict[assembly] + chromosome_str + ".fa"
        if assembly_zipped:
            url = url + ".gz"
        if not os.path.exists(local_assembly_dir):
            os.makedirs(local_assembly_dir)
        wget.download(url, out=local_assembly_dir)

        downloaded_filepath = local_assembly_dir + "/" + url.rsplit('/', 1)[-1]
        if assembly_zipped:
            with gzip.open(downloaded_filepath, 'rb') as f_in:
                with open(local_filepath, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        elif removed_chr == True:
            os.rename(downloaded_filepath, local_filepath)

    def get(self, start, end):
        start = start - 1 
        return self.fa[self.chromosome][start:end].seq
