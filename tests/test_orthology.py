# Paulo Nuin March 2020

import app
from click.testing import CliRunner

import sys
import os
import glob
import requests
import jsonschema
import simplejson as json
import pytest

sys.path.append('../src')

JSON_FILES = '../json'
OUTPUT_DIR = '../output'
SCHEMAS = '../schemas'
RELEASE_VERSION = '3.0.0'


def test_generate_json_file():
    """

    :return:
    """

    print('Starting test')
    runner = CliRunner()
    result = runner.invoke(app.main, ['--orthology'])
    print(result)
    print('File generated')

    assert result.exit_code == 0


def test_validate_generated_file():
    """

    :return:
    """

    print('Starting test')
    print('Checking file ' + OUTPUT_DIR + '/' + 'agr_orthologs-' + RELEASE_VERSION + '.json')

    if not os.path.isfile(OUTPUT_DIR + '/' + 'agr_orthologs-' + RELEASE_VERSION + '.json'):
        assert False
    else:
        with open(SCHEMAS + '/orthology.schema', 'r') as f:
            schema_data = f.read()
        schema = json.loads(schema_data)

        with open(OUTPUT_DIR + '/' + 'agr_orthologs-' + RELEASE_VERSION + '.json') as f:
            result_data = f.read()
        result_file = json.loads(result_data)

        assert jsonschema.validate(result_file, schema) is None


@pytest.mark.skip(reason="not a current test, needs work")
def test_get_json_file():
    """

    :return:
    """

    url = 'http://download.alliancegenome.org/3.0.0/ORTHOLOGY-ALLIANCE-JSON/COMBINED/ORTHOLOGY-ALLIANCE-JSON_COMBINED_17.json'
    local_file = 'ORTHOLOGY-ALLIANCE-JSON_COMBINED_17.json'
    print(glob.glob(JSON_FILES + '/' + local_file))
    if len(glob.glob(JSON_FILES + '/' + local_file)) == 1:
        print('File exists')
    else:
        with open(JSON_FILES + '/' + local_file, 'wb') as f:
            response = requests.get(url)
            f.write(response.content)
#
#     # TODO add assert for file existing or generated/ for files downloaded
#     # assert
