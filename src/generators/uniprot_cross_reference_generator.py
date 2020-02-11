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

    def _read_json_file(self, json_file_location):
        logger.info('Reading JSON file')
        json_file = json.load(json_file_location)
        logger.info('File read')
        return json_file

    def _search_relationships(self, json_file):
        relationships = {}
        for item in json_file:
            if item['GlobalCrossReferenceID'].find('UniProt') >= 0:
                relationships[item['GlobalCrossReferenceID'].replace('UniProtKB:', '')] = item['GeneID']
        return relationships

    def _write_uniprot_file(self, relationships, tab_file):
        output_file = open(tab_file, 'w')
        output_file.write(cls.file_header + '\n')
        output_file.close()

        with open(tab_file, 'a') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            for item in relationships:
                writer.writerow(item, relationships[item])


    def generate_file(self, upload_flag=False):
        logger.info('I will genarate the files')