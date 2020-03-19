import os

import yaml
import logging
import subprocess


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


def compress(cmd):

    logger.info('Running ' + cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()
    stdout, stderr = process.communicate()


    return stdout, stderr, process.returncode