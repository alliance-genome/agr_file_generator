import os

import requests
import yaml
import logging
import subprocess
from collections import OrderedDict

from data_source import DataSource

logger = logging.getLogger(__name__)


# Common configuration variables used throughout the script.
class ContextInfo(object):

    def __init__(self):
        if os.path.exists('config.yaml'):
            config_file = open('config.yaml', 'r')
        else:
            config_file = open('src/config.yaml', 'r')

        self.config = yaml.load(config_file, Loader=yaml.FullLoader)

        # Look for ENV variables to replace default variables from config file.
        for key in self.config.keys():
            try:
                self.config[key] = os.environ[key]
            except KeyError:
                logger.info('Environmental variable not found for \'{}\'. Using config.yaml value.'.format(key))
                pass  # If we don't find an ENV variable, keep the value from the config file.

        logger.debug('Initialized with config values: {}'.format(self.config))


def run_command(cmd):
    logger.info('Running ' + cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()
    stdout, stderr = process.communicate()

    return stdout, stderr, process.returncode


def get_neo_uri(config_info):
    if config_info.config['NEO4J_HOST']:
        uri = "bolt://" + config_info.config['NEO4J_HOST'] + ":" + str(config_info.config['NEO4J_PORT'])
        if config_info.config['DEBUG']:
            logger.info("Using db URI: {}".format(uri))
        return uri
    else:
        logger.error("Need to set NEO4J_HOST env variable")
        exit()


def get_ordered_species_dict(config_info, taxon_ids):
    species_query = """MATCH (s:Species)
                       RETURN s
                       ORDER BY s.phylogeneticOrder"""
    species_data_source = DataSource(get_neo_uri(config_info), species_query)
    species = OrderedDict()
    for record in species_data_source:
        if record["s"]["primaryKey"] in taxon_ids:
            species[record["s"]["primaryKey"]] = record["s"]["name"]

    return species


def ordered_taxon_species_map_from_data_dictionary(taxon_ids):
    species_url = 'https://raw.githubusercontent.com/alliance-genome/agr_schemas/master/ingest/species/species.yaml'
    logger.info('Reading in ' + species_url)
    response = requests.get(species_url)

    if response.status_code == 200:
        species_yaml = yaml.load(response.content, Loader=yaml.FullLoader)
        species = OrderedDict()
        for species_obj in sorted(species_yaml, key=lambda x: x['phylogenicOrder']):
            if species_obj['taxonId'] in taxon_ids:
                species[species_obj['taxonId']] = species_obj['fullName']

        return species
    else:
        logger.critical('unable to download ' + species_url + ' with status code: ' + str(response.status_code))
        exit(-1)
