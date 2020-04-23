"""
.. module:: gene_cross_reference_file_generators
    :platform: any
    :synopsis: Module that generates the Gene Cross Reference files
.. moduleauthor:: AGR consrotium

"""

import os
import logging
import csv
import json
from time import gmtime, strftime

from upload import upload
from .header import create_header

logger = logging.getLogger(name=__name__)


class GeneCrossReferenceFileGenerator:
    """
    TBA
    """

    def __init__(self, gene_cross_references, generated_files_folder, config_info):
        """

        :param gene_cross_references:
        :param generated_files_folder:
        :param config_info:
        """
        self.gene_cross_references = gene_cross_references
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _generate_header(cls, config_info):
        """

        :param config_info:
        :return:
        """

        return create_header('Gene Cross Reference', config_info.config['RELEASE_VERSION'],
                             taxon_ids="# TaxonIDs: " + taxon_ids,
                             species='')


    def generate_file(self, upload_flag=False):
        """

        :param upload_flag:
        :return:
        """
        TSVfilename = 'agr-gene-cross-references-' + self.config_info.config['RELEASE_VERSION'] + '.tsv'
        JSONfilename = 'agr-gene-cross-references-json-' + self.config_info.config['RELEASE_VERSION'] + '.json'
        output_filepath = os.path.join(self.generated_files_folder, TSVfilename)
        output_filepath_json = os.path.join(self.generated_files_folder, JSONfilename)
        gene_cross_reference_file = open(output_filepath, 'w')
        gene_cross_reference_file.write(self._generate_header(self.config_info))
        columns = ['GeneID',
                   'GlobalCrossReferenceID',
                   'CrossReferenceCompleteURL',
                   'ResourceDescriptorPage',
                   'TaxonID']

        tsv_writer = csv.DictWriter(gene_cross_reference_file, delimiter='\t', fieldnames=columns, lineterminator="\n")
        tsv_writer.writeheader()
        listofxrefs = []
        for data in self.gene_cross_references:
            listofxrefs.append(data)
            row = dict(zip(columns, [None] * len(columns)))
            row['GeneID'] = data['GeneID']
            row['GlobalCrossReferenceID'] = data['GlobalCrossReferenceID']
            row['CrossReferenceCompleteURL'] = data['CrossReferenceCompleteURL']
            row['ResourceDescriptorPage'] = data['ResourceDescriptorPage']
            row['TaxonID'] = data['TaxonID']
            tsv_writer.writerows([row])

        gene_cross_reference_file.close()

        with open(output_filepath_json, 'w') as outfile:
            json.dump(listofxrefs, outfile)
        outfile.close()

        if upload_flag:
            logger.info("Submitting to FMS")
            process_name = "1"
            logger.info("uploading TSV version of the gene cross references file.")
            upload.upload_process(process_name, TSVfilename, self.generated_files_folder, 'GENECROSSREFERENCE',
                                  'COMBINED', self.config_info)
            logger.info("uploading JSON version of the gene cross references file.")
            upload.upload_process(process_name, JSONfilename, self.generated_files_folder, 'GENECROSSREFERENCEJSON',
                                  'COMBINED', self.config_info)
