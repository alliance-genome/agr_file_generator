from collections import OrderedDict
from itertools import groupby
from operator import itemgetter
import atexit
import logging
import os
import re
import shutil
import tempfile

import pytest

import agr.app as agr_app


logger = logging.getLogger(name=__package__)

VCF_DATA = OrderedDict()

_temp_folders = set()

_line_split_regex = re.compile(r'\W+')


@atexit.register
def cleanup_temp_folders():
    for fldr in  _temp_folders:
        shutil.rmtree(fldr)

# Mapping from generated file name to a set of example data that should
# appear in the generated file for each mod.
EXAMPLE_CASES = {
    'GRCm38-2.0.0.vcf': [{'CHROMO': '13',
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


def parse_generated_file(path, assembly):
    return parse_vcf_file(path)


def make_gen_files_fixture(asm_cached=False):
    global VCF_DATA
    tempdir = tempfile.mkdtemp()
    _temp_folders.add(tempdir)
    gf_folder_path = os.path.join(tempdir, 'generated_files')
    if asm_cached:
        fasta_sequences_folder = 'sequences'
    else:
        fasta_sequences_folder = os.path.join(tempdir, 'sequences')
    logger.debug('FOLDER PATH:', gf_folder_path)
    os.makedirs(gf_folder_path)
    agr_app.main(generated_files_folder=gf_folder_path,
                 fasta_sequences_folder=fasta_sequences_folder)
    for gf in os.listdir(gf_folder_path):
        logger.debug(gf)
        path = os.path.join(gf_folder_path, gf)
        VCF_DATA[path] = parse_vcf_file(path)
    return gf_folder_path


@pytest.fixture(scope='module')
def run_generate_files():
    return make_gen_files_fixture(asm_cached=True)


def check_files_generated(fixture):
    assert VCF_DATA, 'VCF files not generated'
    for (path, records) in VCF_DATA.items():
        assert os.path.isfile(path)
        assert path.endswith('.vcf')


def test_files_generated_asm_cached(run_generate_files):
    """Run the code to genreate files, assumes the assembly sequence FASTA files
    have been pre-downloaded.

    If they have not, this test should still pass but the test will take a lot longer
    to run.
    """
    check_files_generated(run_generate_files)


def test_ids_unique_in_files(run_generate_files):
    for (path, records) in VCF_DATA.items():
        ids = list(rec['ID'] for rec in records)
        assert ids == list(ids), 'Duplicate id in file:' + path


def vcf_data_by_filename_and_id():
    org_data = OrderedDict()
    for (path, records) in VCF_DATA.items():
        data = org_data[path] = OrderedDict()
        for record in records:
            data.setdefault(record['ID'], []).append(record)
    return org_data


def test_example_expectations(run_generate_files):
    vcf_data = vcf_data_by_filename_and_id()
    for path in vcf_data:
        filename = os.path.basename(path)
        examples = EXAMPLE_CASES.get(filename)
        if not examples:
            logger.info('No examples for ', filename, ', skipping..')
            continue
        parsed_vcf = vcf_data.get(filename)
        if parsed_vcf is None:
            continue
        for example in examples:
            vcf_record = parsed_vcf.get(examples['ID'])
            assert vcf_record, 'No matching VCF data found for example with id: ' + examples['ID']
            assert vcf_record == example, 'Mismatch between example and parsed VCF record'


def test_generated_files_sorted_by_chr_and_pos(run_generate_files):
    for (path, records) in VCF_DATA.items():
        select_chromosome = itemgetter('CHROM')
        chromosomes = list(map(select_chromosome, records))
        assert chromosomes == sorted(chromosomes), 'Chromosomes not alphabetically sorted'
        for (chromo, recs) in groupby(records, select_chromosome):
            row = list(recs)
            positions = list(int(col['POS']) for col in row)
            assert positions == sorted(positions), 'Positions are not sorted in correct order'
