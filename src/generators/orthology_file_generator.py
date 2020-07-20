import os
import logging
import json
import csv

import upload
from .header import create_header
from validators import json_validator

logger = logging.getLogger(name=__name__)


class OrthologyFileGenerator:

    def __init__(self, orthologs, generated_files_folder, config_info):
        self.orthologs = orthologs
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _generate_header(cls, config_info, taxon_ids, data_format):

        return create_header('Orthology', config_info.config['RELEASE_VERSION'],
                             stringency_filter='Stringent',
                             config_info=config_info,
                             taxon_ids=taxon_ids,
                             data_format=data_format)

    def generate_file(self, upload_flag=False, validate_flag=False):

        file_basename = "agr_orthologs-" + self.config_info.config['RELEASE_VERSION']
        fields = ["Gene1ID",
                  "Gene1Symbol",
                  "Gene1SpeciesTaxonID",
                  "Gene1SpeciesName",
                  "Gene2ID",
                  "Gene2Symbol",
                  "Gene2SpeciesTaxonID",
                  "Gene2SpeciesName",
                  "Algorithms",
                  "AlgorithmsMatch",
                  "OutOfAlgorithms",
                  "IsBestScore",
                  "IsBestRevScore"]

        processed_orthologs = []
        taxon_ids = set()
        for ortholog in self.orthologs:
            num_algorithms = ortholog["numAlgorithmMatch"] + ortholog["numAlgorithmNotMatched"]
            taxon_ids.add(ortholog["species1TaxonID"])
            taxon_ids.add(ortholog["species2TaxonID"])
            processed_orthologs.append(dict(zip(fields, [ortholog["gene1ID"],
                                                         ortholog["gene1Symbol"],
                                                         ortholog["species1TaxonID"],
                                                         ortholog["species1Name"],
                                                         ortholog["gene2ID"],
                                                         ortholog["gene2Symbol"],
                                                         ortholog["species2TaxonID"],
                                                         ortholog["species2Name"],
                                                         ortholog["Algorithms"],
                                                         str(ortholog["numAlgorithmMatch"]),
                                                         num_algorithms,
                                                         ortholog["best"],
                                                         ortholog["bestRev"]])))

        json_filename = file_basename + ".json"
        json_filepath = os.path.join(self.generated_files_folder, json_filename)
        with open(json_filepath, 'w') as outfile:
            contents = {'metadata': self._generate_header(self.config_info, taxon_ids, 'json'),
                        'data': processed_orthologs}
            json.dump(contents, outfile)

        for processed_ortholog in processed_orthologs:
            processed_ortholog['Algorithms'] = "|".join(set(processed_ortholog['Algorithms']))

        tsv_filename = file_basename + ".tsv"
        tsv_filepath = os.path.join(self.generated_files_folder, tsv_filename)
        tsv_file = open(tsv_filepath, 'w')
        tsv_file.write(self._generate_header(self.config_info, taxon_ids, 'tsv'))

        tsv_writer = csv.DictWriter(tsv_file, delimiter='\t', fieldnames=fields, lineterminator="\n")
        tsv_writer.writeheader()
        tsv_writer.writerows(processed_orthologs)
        tsv_file.close()

        if validate_flag:
            json_validator.JsonValidator(json_filepath, 'orthology').validateJSON()
            if upload_flag:
                logger.info("Submitting orthology filse to FMS")
                process_name = "1"
                upload.upload_process(process_name, json_filepath, self.generated_files_folder, 'ORTHOLOGY-ALLIANCE-JSON', 'COMBINED', self.config_info)
                upload.upload_process(process_name, tsv_filepath, self.generated_files_folder, 'ORTHOLOGY-ALLIANCE', 'COMBINED', self.config_info)
