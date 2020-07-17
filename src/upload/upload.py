# Functions for use in downloading files.
import os

import logging
import requests
# from requests_toolbelt.utils import dump
from retry import retry

logger = logging.getLogger(__name__)


def upload_file(worker, filename, save_path, upload_file_prefix, config_info):
    with open(os.path.join(save_path, filename), 'rb') as fp:
        file_to_upload = {upload_file_prefix: fp}
        logger.info(file_to_upload)
        if config_info.config['API_KEY']:
            headers = {
                'Authorization': 'Bearer {}'.format(config_info.config['API_KEY'])
            }
        else:
            headers = {}
        logger.debug('{}: Attempting upload of data file: {}'.format(worker, os.path.join(save_path, filename)))
        logger.debug('{}: Attempting upload with header: {}'.format(worker, headers))
        logger.info("{}: Uploading data to {}) ...".format(worker, config_info.config['FMS_API_URL'] + '/api/data/submit/'))
        try:
            response = requests.post(config_info.config['FMS_API_URL'] + '/api/data/submit', files=file_to_upload, headers=headers)
            logger.info(response.text)
        except requests.exceptions.RequestException as e:
            raise(e)


@retry(tries=5, delay=5, logger=logger)
def upload_process(worker, filename, save_path, data_type, data_sub_type, config_info):

    upload_file_prefix = '{}_{}_{}'.format(config_info.config['RELEASE_VERSION'], data_type, data_sub_type)

    # Attempt to grab MD5 for the latest version of the file.
    logger.info(config_info.config['FMS_API_URL'] + '/api/datafile/by/{}/{}?latest=true'.format(data_type, data_sub_type))
    upload_file(worker, filename, save_path, upload_file_prefix, config_info)
