import ntpath
import logging
import collections
from operator import itemgetter
from collections import OrderedDict
from itertools import groupby


logger = logging.getLogger(name=__name__)


EXAMPLE_CASES = {
  'GRCm38': [{'CHROMO': '13',
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
              'INFO': ''}]}


class VcfValidator:

    def __init__(self, filepath):
        self.filepath = filepath
        self.filename = ntpath.basename(filepath)
        self.assembly = self.filename.split('-')[0]

    def parse_vcf_file(self, path):
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

    def check_examples(self, parsed_vcf):
        examples = EXAMPLE_CASES.get(self.assembly)
        if not examples:
            logger.info('No examples for ' + self.filename + ', skipping ...')
        else:
            for example in examples:
                vcf_record = parsed_vcf.get(examples['ID'])
                if vcf_record not in examples['ID']:
                    logger.error('No matching VCF data found for example with id: ' + examples['ID'])
                    exit(-1)
                if vcf_record != example:
                    logger.error('Mismatch between example and parsed VCF record')

    def check_sorted_by_chromosome_and_position(self, parsed_vcf):
        select_chromosome = itemgetter('CHROM')
        chromosomes = list(map(select_chromosome, parsed_vcf))
        assert chromosomes == sorted(chromosomes), 'Chromosomes not alphabetically sorted'
        for (chromo, recs) in groupby(parsed_vcf, select_chromosome):
            row = list(recs)
            positions = list(int(col['POS']) for col in row)
            if positions != sorted(positions):
                logger.error('Positions are not sorted in correct order')
                exit(-1)

        logger.info('Sorted by chromosome and position')

    def check_duplicate_entries(self, parsed_vcf):
        all_entries = []
        for variant in parsed_vcf:
            all_entries.append(variant['ID'])

        print(all_entries)
        if len(all_entries) == len(set(all_entries)):
            logger.info("No duplicate enteries")
        else:
            logger.error("At least one Duplicate entery")
            logger.error([item for item, count in collections.Counter(all_entries).items() if count > 1])
            exit(-1)

    def validate_vcf(self):
        logger.info("Validating VCF: %s" % self.filename)

        parsed_vcf = self.parse_vcf_file(self.filepath)

        self.check_examples(parsed_vcf)
        self.check_sorted_by_chromosome_and_position(parsed_vcf)
        self.check_duplicate_entries(parsed_vcf)
