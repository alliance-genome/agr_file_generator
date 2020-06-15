import app
from click.testing import CliRunner
import logging
import os
import glob
import sys
from collections import OrderedDict
from operator import itemgetter
from itertools import groupby

sys.path.append('../src')

logger = logging.getLogger(name=__package__)

OUTPUT_DIR = '../output'
RELEASE_VERSION = '3.0.0'
VCF_DATA = OrderedDict()

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


def parse_vcf_file(path):
    """

    :param path:
    :return:
    """

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


def vcf_data_by_filename_and_id():
    """

    :return:
    """

    VCF_DATA = get_full_vcf_data()
    org_data = OrderedDict()

    for (path, records) in VCF_DATA.items():
        print(path)
        data = org_data[path] = OrderedDict()
        for record in records:
            data.setdefault(record['ID'], []).append(record)

    return org_data


def get_full_vcf_data():
    """

    :return:
    """

    vcf_data = OrderedDict()
    for gf in glob.glob(OUTPUT_DIR + '/*.vcf'):
        path = os.path.join(gf)
        vcf_data[path] = parse_vcf_file(path)

    return vcf_data


def test_file_generation():
    """

    :return:
    """

    runner = CliRunner()
    result = runner.invoke(app.main, ['--vcf'])
    assert result.exit_code == 0


def test_generated_files():
    """

    :return:
    """

    get_full_vcf_data()
    for (path, records) in VCF_DATA.items():
        print(path)
        assert os.path.isfile(path)
        assert path.endswith('.vcf')


def test_example_expectations():
    """

    :return:
    """

    vcf_data = vcf_data_by_filename_and_id()
    print(len(vcf_data))
    for path in vcf_data:
        filename = os.path.basename(path)
        print(filename)
        examples = EXAMPLE_CASES.get(filename)
        if not examples:
            print('No examples for ', filename, ', skipping ...')
            continue
        parsed_vcf = vcf_data.get(filename)
        if parsed_vcf is None:
            continue
        for example in examples:
            vcf_record = parsed_vcf.get(examples['ID'])
            assert vcf_record, 'No matching VCF data found for example with id: ' + examples['ID']
            assert vcf_record == example, 'Mismatch between example and parsed VCF record'


def test_generated_files_sorted_by_chr_and_pos():
    """

    :return:
    """

    vcf_data = get_full_vcf_data()
    for (path, records) in vcf_data.items():
        select_chromosome = itemgetter('CHROM')
        chromosomes = list(map(select_chromosome, records))
        assert chromosomes == sorted(chromosomes), 'Chromosomes not alphabetically sorted'
        for (chromo, recs) in groupby(records, select_chromosome):
            row = list(recs)
            positions = list(int(col['POS']) for col in row)
            assert positions == sorted(positions), 'Positions are not sorted in correct order'


def test_duplicate_entries():
    """

    :return:
    """

    vcf_data = vcf_data_by_filename_and_id()
    all_entries = []
    for i in vcf_data:
        for key, value in vcf_data[i].items():
            all_entries.append(key)
    assert len(all_entries) == len(set(all_entries))
