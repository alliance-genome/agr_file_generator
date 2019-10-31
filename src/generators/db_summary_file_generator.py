import os
import time
import json
import logging

import upload


logger = logging.getLogger(name=__name__)


class DbSummaryFileGenerator:

    empty_value_marker = '.'

    def __init__(self, entities, generated_files_folder, config_info):
        self.entitites = entities
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder


    def generate_file(self, upload_flag=False):
        entities = dict()
        for record in self.entitites:

            frequency = record["frequency"]
            entityTypes = record["entityTypes"]
            if (len(entityTypes) == 1 and entityTypes[0] != "Load"):
                entities[entityTypes[0]] = frequency
            elif len(entityTypes) == 2:
                if entityTypes[1] in entities:
                    entities[entityTypes[1]][entityTypes[0]] = frequency
                else:
                    entities[entityTypes[1]] = {entityTypes[0]: frequency}
                if (entityTypes[0] not in entities) or (entityTypes[0] in entities and isinstance(entities[entityTypes[0]], dict)):
                    if entityTypes[0] in entities:
                        entities[entityTypes[0]][entityTypes[1]] = frequency
                    else:
                        entities[entityTypes[0]] = {entityTypes[1]: frequency}

         

        entityKeys = list(entities.keys()).copy()
        for key in entityKeys:
            if not isinstance(entities[key], int):
                if len(entities[key].keys()) == 1:
                    subKey = list(entities[key].keys())[0]
                    if subKey in entities:
                        del entities[key]
    
        summary = { "overview": entities }

        filename =  'db-summary-' + self.config_info.config['RELEASE_VERSION'] + '.json'
        filepath = os.path.join(self.generated_files_folder, filename)
        print(filepath)
        with open(filepath, 'w') as json_file:
            json.dump(summary, json_file, sort_keys=True, indent=4)

        if upload_flag:
            logger.info("Submitting to FMS")
            process_name = "1"
            upload.upload_process(process_name, filename, self.generated_files_folder, 'DBSUMMARY', config_info.config['RELEASE_VERSION'], self.config_info)
