

import logging
import os
import glob
import sys
from collections import OrderedDict
import pytest
sys.path.append('../src')

from common import ContextInfo
import app
import click
from click.testing import CliRunner

logger = logging.getLogger(name=__package__)

OUTPUT_DIR = '../output'
RELEASE_VERSION = '3.0.0'
VCF_DATA = OrderedDict()

EXAMPLE_FILE = 'GRCm38-' + RELEASE_VERSION

EXAMPLE_CASES = {
    EXAMPLE_FILE: [{'CHROMO': '13',
                          'POS': '50540171',
                          'ID': 'ZFIN:ZDB-ALT-160601-8105',
                         'REF': 'C',
                          'ALT': 'T',
                          'QUAL': '',
                          'FILTER': '',
                          'INFO': ''},
                         {'CHROMO': '5',
                          'POS': '72118556',
                          'ID': 'ZFIN:ZDB-ALT-170321-11',
                          'REF': 'CGCTTTGA',
                          'ALT': 'C',
                          'QUAL': '',
                          'FILTER': ''},
                         {'CHROM': '10',
                          'POS': '16027812',
                          'ID': 'ZFIN:ZDB-ALT-180207-16',
                          'REF': 'G',
                          'ALT': 'GCCGTT',
                          'QUAL': '',
                          'FILTER': '',
                          'INFO': ''}]
}

def parse_vcf_file(path):
    data = []
    headers = []
    with open(path, 'rt') as fp:
        for line in fp:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                cols = line[1:].split('\t')
                headers.extend(cols)
                continue
            else:
                cols = line.split('\t')
            data.append(OrderedDict(zip(headers, cols)))
    return data


def test_file_generation():

    runner = CliRunner()
    result = runner.invoke(app.main, ['--vcf'])
    assert result.exit_code == 0


def test_generated_files():

    for gf in glob.glob(OUTPUT_DIR + '/*.vcf'):
        path = os.path.join(gf)
        VCF_DATA[path] = parse_vcf_file(path)

    for (path, records) in VCF_DATA.items():
        print(path)
        assert os.path.isfile(path)
        assert path.endswith('.vcf')










    # runner = CliRunner()
    # result = runner.invoke(app.main, ['--vcf'])
    # print(result)
    # print('File generated ... or not')
    #
    # assert result.exit_code == 0
    #
    # assert VCF_DATA, 'VCF files not generated. Please generate files and re-run tests'
    # for (path, records) in VCF_DATA.items():
    #     assert os.path.isfile(path)
    #     assert path.endswith('.vcf')


