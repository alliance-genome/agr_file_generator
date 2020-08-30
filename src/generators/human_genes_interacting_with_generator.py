import os
import logging
import json
import csv

import upload
from .header import create_header
from validators import json_validator

logger = logging.getLogger(name=__name__)


class HumanGenesInteractingWithGenerator:

    def __init__(self, interactions, config_info, generated_files_folder):
        self.interactions = interactions
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _generate_header(cls, config_info, data_format):

        return create_header('Human genes encoding proteins that interact with SARS-CoV-2 proteins', config_info.config['RELEASE_VERSION'],
                             config_info=config_info,
                             taxon_ids=["NCBI:txid9606"],
                             data_format=data_format)

    def generate_file(self, upload_flag=False, validate_flag=False):
        file_basename = "agr-human_genes_interacting_with-" + self.config_info.config['RELEASE_VERSION']
        fields = ["GeneID",
                  "GeneSymbol",
                  "Name"]

        processed_interactions = []
        for interaction in self.interactions:
            processed_interactions.append(dict(zip(fields, [interaction["GeneID"],
                                                            interaction["Symbol"],
                                                            interaction["Name"]])))

        json_filename = file_basename + ".json"
        json_filepath = os.path.join(self.generated_files_folder, json_filename)
        with open(json_filepath, 'w') as outfile:
            contents = {'metadata': self._generate_header(self.config_info, 'json'),
                        'data': processed_interactions}
            json.dump(contents, outfile)

        tsv_filename = file_basename + ".tsv"
        tsv_filepath = os.path.join(self.generated_files_folder, tsv_filename)
        tsv_file = open(tsv_filepath, 'w')
        tsv_file.write(self._generate_header(self.config_info, 'tsv'))

        tsv_writer = csv.DictWriter(tsv_file, delimiter='\t', fieldnames=fields, lineterminator="\n")
        tsv_writer.writeheader()
        tsv_writer.writerows(processed_interactions)
        tsv_file.close()

        if validate_flag:
            json_validator.JsonValidator(json_filepath, 'human-genes-interacting-with').validateJSON()
            if upload_flag:
                logger.info("Submitting human genes interacting with filse to FMS")
                process_name = "1"
                upload.upload_process(process_name, json_filepath, self.generated_files_folder, 'Human-genes-interacting-with-JSON', 'SARS-CoV-2', self.config_info)
                upload.upload_process(process_name, tsv_filepath, self.generated_files_folder, 'Human-genes-interacting-with', 'SARS-CoV-2', self.config_info)
