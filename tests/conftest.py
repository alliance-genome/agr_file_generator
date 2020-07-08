# Paulo Nuin June 2020

# tesr configuration will be added here

from collections import OrderedDict
import glob
import os
import pytest

OUTPUT_DIR = '../output'



def test_connection():

    pass


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

@pytest.fixture()
def get_full_vcf_data():
    """

    :return:
    """

    vcf_data = OrderedDict()
    for gf in glob.glob(OUTPUT_DIR + '/*.vcf'):
        path = os.path.join(gf)
        vcf_data[path] = parse_vcf_file(path)

    return vcf_data


@pytest.fixture()
def vcf_data_by_filename_and_id(get_full_vcf_data):
    """

    :return:
    """

    # VCF_DATA = get_full_vcf_data()
    org_data = OrderedDict()

    for (path, records) in get_full_vcf_data.items():
        print(path)
        data = org_data[path] = OrderedDict()
        for record in records:
            data.setdefault(record['ID'], []).append(record)

    return org_data