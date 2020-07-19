"""
.. module:: db_summary_file_generator
    :platform: any
    :synopsis: Module that generates a summary of the DB entities
.. moduleauthor:: AGR consortium

"""

import os
# import time
import json
import logging

import upload


logger = logging.getLogger(name=__name__)


class DbSummaryFileGenerator:
    """
    TBA
    """

    empty_value_marker = '.'

    def __init__(self, entities, generated_files_folder, config_info):
        """

        :param entities:
        :param generated_files_folder:
        :param config_info:
        """
        self.entitites = entities
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder

    def __get_entity_counters(self):
        """Get the count for each entity in database"""

        entity_counters = dict()
        for record in self.entitites:

            frequency = record["frequency"]
            entity_types = record["entityTypes"]
            if (len(entity_types) == 1 and entity_types[0] != "Load"):
                entity_counters[entity_types[0]] = frequency
            elif len(entity_types) == 2:
                if entity_types[1] in entity_counters and isinstance(entity_types[1], dict):
                    entity_counters[entity_types[1]][entity_types[0]] = frequency
                else:
                    entity_counters[entity_types[1]] = {entity_types[0]: frequency}
                if (entity_types[0] not in entity_counters) or (entity_types[0] in entity_counters and isinstance(entity_counters[entity_types[0]], dict)):
                    if entity_types[0] in entity_counters:
                        entity_counters[entity_types[0]][entity_types[1]] = frequency
                    else:
                        entity_counters[entity_types[0]] = {entity_types[1]: frequency}

        for key in list(entity_counters.keys()).copy():
            if not isinstance(entity_counters[key], int):
                if len(entity_counters[key].keys()) == 1:
                    subKey = list(entity_counters[key].keys())[0]
                    if subKey in entity_counters:
                        del entity_counters[key]

        return entity_counters

    def __generate_overview(self):
        """Generate overview counts for each nod label and sub label"""

        entity_counters = self.__get_entity_counters()

        overview = []
        for node_label in list(entity_counters.keys()).copy():
            if isinstance(entity_counters[node_label], dict):
                counter = 0
                for sub_label, sub_node_counter in entity_counters[node_label].items():
                    overview.append({"nodeLabel": sub_label,
                                     "parentLabel": node_label,
                                     "counter": sub_node_counter})
                    counter += sub_node_counter

                overview.append({"nodeLabel": node_label,
                                 "counter": counter})
            else:
                overview.append({"nodeLabel": node_label,
                                 "counter": entity_counters[node_label]})

        return overview

    def generate_file(self, upload_flag=False, validate_flag=False):
        """

        :param upload_flag:
        :return:
        """

        summary = {"overview": self.__generate_overview()}

        filename = 'db-summary-' + self.config_info.config['RELEASE_VERSION'] + '.json'
        filepath = os.path.join(self.generated_files_folder, filename)
        logger.info(filepath)
        with open(filepath, 'w') as json_file:
            json.dump(summary, json_file, sort_keys=True, indent=4)

        if validate_flag:
            if upload_flag:
                logger.info("Submitting to FMS")
                process_name = "1"
                upload.upload_process(process_name,
                                      filename,
                                      self.generated_files_folder,
                                      'DB-SUMMARY',
                                      self.config_info.config['RELEASE_VERSION'],
                                      self.config_info)
