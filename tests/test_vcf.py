
from click.testing import CliRunner
import logging
import os
import sys
import pytest
from collections import OrderedDict
from operator import itemgetter
from itertools import groupby

sys.path.append('../src')
import app

logger = logging.getLogger(name=__package__)

OUTPUT_DIR = '../output'
RELEASE_VERSION = '3.1.0'
EXAMPLE_FILE = 'GRCm38-' + RELEASE_VERSION + '.vcf'
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

@pytest.mark.skip()
def test_file_generation():
    """

    :return:
    """

    runner = CliRunner()
    result = runner.invoke(app.main, ['--vcf'])
    assert result.exit_code == 0


def test_generated_files(get_full_vcf_data):
    """

    :return:
    """

    for (path, records) in get_full_vcf_data.items():
        print(path)
        assert os.path.isfile(path)
        assert path.endswith('.vcf')


def test_example_expectations(vcf_data_by_filename_and_id):
    """

    :return:
    """

    print(len(vcf_data_by_filename_and_id))
    for path in vcf_data_by_filename_and_id:
        filename = os.path.basename(path)
        print(filename)
        examples = EXAMPLE_CASES.get(filename)
        if not examples:
            print('No examples for ', filename, ', skipping ...')
            continue
        parsed_vcf = vcf_data_by_filename_and_id.get(filename)
        if parsed_vcf is None:
            continue
        for example in examples:
            vcf_record = parsed_vcf.get(examples['ID'])
            assert vcf_record, 'No matching VCF data found for example with id: ' + examples['ID']
            assert vcf_record == example, 'Mismatch between example and parsed VCF record'


def test_generated_files_sorted_by_chr_and_pos(get_full_vcf_data):
    """

    :return:
    """

    for (path, records) in get_full_vcf_data.items():
        select_chromosome = itemgetter('CHROM')
        chromosomes = list(map(select_chromosome, records))
        assert chromosomes == sorted(chromosomes), 'Chromosomes not alphabetically sorted'
        for (chromo, recs) in groupby(records, select_chromosome):
            row = list(recs)
            positions = list(int(col['POS']) for col in row)
            assert positions == sorted(positions), 'Positions are not sorted in correct order'


def test_duplicate_entries(vcf_data_by_filename_and_id):
    """

    :return:
    """


    all_entries = []
    for i in vcf_data_by_filename_and_id:
        for key, value in vcf_data_by_filename_and_id[i].items():
            all_entries.append(key)
    assert len(all_entries) == len(set(all_entries))
