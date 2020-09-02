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
        self.assembly = assembly
        self.alleles = alleles
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder

    def _get_vcf_start_position(variant):
        so_term = variant['soTerm']

        if so_term == 'deletion':
            return variant['start'] - 1
        elif so_term in ['insertion', 'point_mutation', 'MNV']:
            return variant['start']
        elif so_term == 'delins':
            if variant['genomicVariantSequence'] == '':
                return variant['start'] - 1
            elif len(variant['genomicVariantSequence']) == len(variant['genomicReferenceSequence']):
                return variant['start']
            else:
                return variant['start'] - 1
        else:
            logger.fatal('New SoTerm that We need to add logic for: %r', so_term)
            return None

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
                variant_rows = []
                allele_start = -1
                allele_end = 0

                for variant in allele["variants"]:
                    start = AlleleGffFileGenerator._get_vcf_start_position(variant)
                    end = variant['end']
                    if start < allele_start or allele_start == -1:
                        allele_start = start
                    if end > allele_end:
                        allele_end = end

                    gene_ids = []
                    gene_symbols = []
                    gene_impacts = []
                    gene_consequences = []
                    for glc in variant['geneLevelConsequences']:
                        gene_ids.append(glc['geneID'])
                        gene_symbols.append(glc['geneSymbol'])
                        gene_impacts.append(glc['impact'])
                        gene_consequences.append(glc['geneLevelConsequence'].replace(",", "|"))

                    column_nine = ';'.join(['Parent=' + allele['ID'],
                                            'geneID=' + ','.join(gene_ids),
                                            'geneSymbol=' + ','.join(gene_symbols),
                                            'geneImpacts=' + ','.join(gene_impacts),
                                            'geneConsequences=' + ','.join(gene_consequences)])

                    variant_rows.append([variant['chromosome'],
                                         '.',
                                         variant['soTerm'],
                                         str(start),
                                         str(end),
                                         '.', '.', '.',
                                         column_nine])

                column_nine = ';'.join(['ID=' + allele['ID'],
                                        'symbol=' + allele['symbol'],
                                        'symbol_text=' + allele['symbol_text']])
                allele_file.write("\t".join([allele['chromosome'],
                                             '.',
                                             'biological_region',
                                             str(allele_start),
                                             str(allele_end),
                                             '.', '.', '.',
                                             column_nine]) + "\n")
                for row in variant_rows:
                    allele_file.write('\t'.join(row) + '\n')

        if validate_flag:
            process_name = "1"
            if upload_flag:
                logger.info("Submitting Allele GFF (" + self.assembly + ") to FMS")
                upload.upload_process(process_name,
                                      filename,
                                      self.generated_files_folder,
                                      'ALLELE-GFF',
                                      self.assembly.replace('_', '').replace('.', ''),
                                      self.config_info)
