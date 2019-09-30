import os
import logging
from time import gmtime, strftime

#from upload.upload import upload_process
import upload

logger = logging.getLogger(name=__name__)


class OrthologyFileGenerator:

    file_header_template = """#########################################################################
#
# Ortholog File
# Source: Alliance of Genome Resources (Alliance)
# Filter: stringent
# Filter Details:
#   Ortholog needs to be called by at least 3 Algorithms
#        or is being called by either ZFIN or HGNC algorithms
#        and either is best score or is best reverse score
#   or ortholog is called by 2 algorithms and is best score and best reverse score
# Datebase Version: {databaseVersion}
# Date: {datetimeNow}
#
#########################################################################
"""

    def __init__(self, orthologs, generated_files_folder, config_info):
        self.orthologs = orthologs
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _generate_header(cls, config_info):
        return cls.file_header_template.format(datetimeNow=strftime("%Y-%m-%d %H:%M:%S", gmtime()),
                                               databaseVersion=config_info.config['RELEASE_VERSION'])

    def generate_file(self, upload_flag=False):
        filename = "agr_orthologs-" + self.config_info.config['RELEASE_VERSION'] + ".tsv"
        outputFilepath = os.path.join(self.generated_files_folder, filename)
        orthology_file = open(outputFilepath,'w')
        orthology_file.write(self._generate_header(self.config_info))
        columns = ["Gene1ID",
                   "Gene1Symbol",
                   "Gene1SpeciesTaxonID",
                   "Gene1SpeciesName",
                   "Gene2ID",
                   "Gene2Symbol",
                   "Gene1SpeciesTaxonID",
                   "Gene1SpeciesName",
                   "Algorithms",
                   "AlgorithmsMatch",
                   "OutOfAlgorithms",
                   "IsBestScore",
                   "IsBestRevScore"]
        orthology_file.write("\t".join(columns) + "\n")

        for ortholog in self.orthologs:
            algorithms = "|".join(set(ortholog["Algorithms"]))
            numAlgorithm = ortholog["numAlgorithmMatch"] + ortholog["numAlgorithmNotMatched"]
            orthology_file.write("\t".join([ortholog["gene1ID"],
                                            ortholog["gene1Symbol"],
                                            ortholog["species1TaxonID"],
                                            ortholog["species1Name"],
                                            ortholog["gene2ID"],
                                            ortholog["gene2Symbol"],
                                            ortholog["species2TaxonID"],
                                            ortholog["species2Name"],
                                            algorithms,
                                            str(ortholog["numAlgorithmMatch"]),
                                            str(numAlgorithm),
                                            ortholog["best"],
                                            ortholog["bestRev"]]) + "\n")
        orthology_file.close()
        if upload_flag:
            logger.info("Submitting to FMS")
            process_name = "1"
            upload.upload_process(process_name, filename, self.generated_files_folder, 'ORTHOLOGY', 'ALL', self.config_info)
