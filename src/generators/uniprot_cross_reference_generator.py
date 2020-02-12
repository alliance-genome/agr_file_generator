import os
import time

import sys
import json
import csv
import logging
import coloredlogs


logger = logging.getLogger(name=__name__)


class UniProtGenerator:

    file_header = '''UniProtID\tIdentifier'''

    def __init__(self, config_info, generated_files_folder, input_folder):
        # self.relationships = relationships
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder
        self.input_folder = input_folder

    def _read_json_file(cls, json_file_location, json_filename):
        logger.info('Reading JSON file')
        json_file = json.load(open(json_file_location + '/' + json_filename))
        logger.info('File read')
        return json_file

    def _search_relationships(cls, json_data):
        logger.info('Looking for relationships')
        relationships = {}
        for item in json_data:
            if item['GlobalCrossReferenceID'].find('UniProt') >= 0:
                relationships[item['GlobalCrossReferenceID'].replace('UniProtKB:', '')] = item['GeneID']
        return relationships

    def _write_uniprot_file(cls, relationships, tab_file):
        logger.info('Writing output file')
        output_file = open(cls.generated_files_folder + '/' + tab_file, 'w')
        output_file.write(cls.file_header + '\n')
        output_file.close()

        with open(cls.generated_files_folder + '/' + tab_file, 'a') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            for item in relationships:
                writer.writerow([item, relationships[item]])
        tsvfile.close()

    def generate_file(self, upload_flag=False):
        json_input = self._read_json_file(self.input_folder, 'GENECROSSREFERENCEJSON_COMBINED_4.json')
        relationships = self._search_relationships(json_input)
        self._write_uniprot_file(relationships, 'CROSSREFERENCEUNIPROT_COMBINED.tsv')
        logger.info('File created')