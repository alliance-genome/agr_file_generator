import os
import sys
import logging
import upload
from headers import create_header

sys.path.append('../')

logger = logging.getLogger(name=__name__)


class AlleleGffFileGenerator:

    empty_value_marker = '.'

    def __init__(self, assembly, alleles, generated_files_folder, config_info):
        self.alleles = alleles
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder

    def generate_assembly_file(self, upload_flag=False, validate_flag=False):
        filename = self.assembly + '-' + self.config_info.config['RELEASE_VERSION'] + '.allele.gff'
        filepath = os.path.join(self.generated_files_folder, filename)
        logger.info('Generating Allele GFF File for assembly %r', self.assembly)
        with open(filepath, 'w') as allele_file:
            header = create_header('Allele GFF',
                                   self.config_info.config['RELEASE_VERSION'],
                                   assembly=self.assembly,
                                   config_info=self.config_info,
                                   data_format='GFF')

            allele_file.write(header)
            for allele in self.alleles:
                print(allele)

        if validate_flag:
            process_name = "1"
            if upload_flag:
                logger.info("Submitting Allele GFF (" + self.assembly + ") to FMS")
                upload.upload_process(process_name, filename, self.generated_files_folder, 'ALLELE-GFF', self.assembly, self.config_info)
