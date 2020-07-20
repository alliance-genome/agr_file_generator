import ntpath
import logging
import collections
from operator import itemgetter
from collections import OrderedDict
from itertools import groupby
from common import run_command

logger = logging.getLogger(name=__name__)


EXAMPLE_CASES = {
  'GRCz11': [{'CHROM': '13',
              'POS': '50540171',
              'ID': 'NC_007124.7:g.50540171C>T',
              'REF': 'C',
              'ALT': 'T',
              'QUAL': '.',
              'FILTER': '.'},
             {'CHROM': '5',
              'POS': '72118556',
              'ID': 'NC_007116.7:g.72118557_72118563del',
              'REF': 'CACTTTGA',
              'ALT': 'C',
              'QUAL': '.',
              'FILTER': '.'},
             {'CHROM': '10',
              'POS': '16027812',
              'ID': 'NC_007121.7:g.16027812_16027813insCCGTT',
              'REF': 'G',
              'ALT': 'GCCGTT',
              'QUAL': '.',
              'FILTER': '.'}
             ]}


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
            logger.info('Checking against examples')
            found = False
            for example in examples:
                for vcf_record in parsed_vcf:
                    if vcf_record['ID'] == example['ID']:
                        found = True
                        for key in example.keys():
                            if example[key] != vcf_record[key]:
                                logger.error('Mismatch between example and parsed VCF record')
                                logger.error("key mismatch: " + key)
                                logger.error("example value: " + example[key])
                                logger.error("VCF record value: " + vcf_record[key])
                                exit(-1)
                if not found:
                    logger.error('No matching VCF data found for example with id: ' + example['ID'])
                    exit(-1)

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

        if len(all_entries) == len(set(all_entries)):
            logger.info("No duplicate enteries")
        else:
            logger.error("At least one Duplicate entery")
            logger.error([item for item, count in collections.Counter(all_entries).items() if count > 1])
            exit(-1)

    def run_vcf_validator_cmd(filepath):
        stdout, stderr, return_code = run_command('vcf-validator ' + filepath)
        logger.info(stdout)
        logger.info(stderr.decode("utf-8"))
        if return_code != 0:
            logger.error("vcf_validate caught error in file")
            exit(-1)


    def validate_vcf(self):
        logger.info("Validating VCF: %s" % self.filename)

        parsed_vcf = self.parse_vcf_file(self.filepath)

        self.check_examples(parsed_vcf)
        self.check_sorted_by_chromosome_and_position(parsed_vcf)
        self.check_duplicate_entries(parsed_vcf)
