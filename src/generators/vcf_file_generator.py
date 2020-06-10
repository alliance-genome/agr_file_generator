import os
import time
import sys
from collections import defaultdict, OrderedDict
from functools import partial
from operator import itemgetter
from common import compress
import logging
import upload

sys.path.append('../')

logger = logging.getLogger(name=__name__)


class VcfFileGenerator:

    empty_value_marker = '.'

    file_header = """##contig=<ID=,length=,assembly={assembly},md5=,species="{species}",taxonomy=x>
##fileDate={datetime}
##fileformat=VCFv4.2
##ALT=<ID=R,Description="IUPAC code R = A/G">
##ALT=<ID=Y,Description="IUPAC code Y = C/T">
##ALT=<ID=S,Description="IUPAC code S = C/G">
##ALT=<ID=W,Description="IUPAC code W = A/T">
##ALT=<ID=M,Description="IUPAC code M = A/C">
##ALT=<ID=K,Description="IUPAC code K = T/G">
##ALT=<ID=B,Description="IUPAC code B = C/G/T">
##ALT=<ID=D,Description="IUPAC code D = A/G/T">
##ALT=<ID=H,Description="IUPAC code H = A/C/T">
##ALT=<ID=V,Description="IUPAC code V = A/C/G">
##INFO=<ID=hgvs_nomenclature,Number=1,Type=String,Description="the HGVS name of the allele">
##INFO=<ID=geneLevelConsequence,Number=.,Type=String,Description="VEP consequence of the variant on the Gene">
##INFO=<ID=transcriptLevelConsequence,Number=.,Type=String,Description="VEP consequence of the variant on the Transcript">
##INFO=<ID=geneImpact,Number=1,Type=String,Description="Variant impact scale for Gene">
##INFO=<ID=transcriptImpact,Number=.,Type=String,Description="Variant impact scale for Transcript">
##INFO=<ID=allele_symbols,Number=.,Type=String,Description="The human readable name of the Allele">
##INFO=<ID=soTerm,Number=1,Type=String,Description="The Sequence Ontology term for the variant">
##INFO=<ID=alleles,Number=.,Type=String,Description="The alleles of the variant">
##INFO=<ID=allele_of_gene_ids,Number=1,Type=String,Description="The gene ids that the Allele is located on">
##INFO=<ID=allele_of_gene_symbols,Number=1,Type=String,Description="The gene names that the Allele is located on">
##INFO=<ID=allele_of_transcript_ids,Number=.,Type=String,Description="The gene ids that the Allele is located on">
##INFO=<ID=allele_of_transcript_gff3_ids,Number=.,Type=String,Description="The transcript gff3ID that the Allele is located on">
##INFO=<ID=allele_symbols_text,Number=.,Type=String,Description="Another human readable representation of the allele">
##phasing=partial
##source=AGR VCF File generator"""

    col_headers = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')

    def __init__(self, variants, generated_files_folder, config_info):
        self.variants = variants
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _add_padded_base_to_variant(cls, variant, so_term):
        padded_base = variant['paddingLeft']
        if not ("," in variant['genomicReferenceSequence'] or "," in variant['genomicVariantSequence']):
            variant['genomicReferenceSequence'] = padded_base + variant['genomicReferenceSequence']
            variant['genomicVariantSequence'] = padded_base + variant['genomicVariantSequence']

    @classmethod
    def _write_vcf_header(cls, vcf_file, assembly, species, config_info):
        dt = time.strftime("%Y%m%d", time.gmtime())
        header = cls.file_header.format(datetime=dt,
                                        database_version=config_info.config['RELEASE_VERSION'],
                                        species=species,
                                        assembly=assembly)
        vcf_file.write(header)
        vcf_file.write('\n')
        vcf_file.write('#' + '\t'.join(cls.col_headers))
        vcf_file.write('\n')

    @classmethod
    def _write_tab_delimited_header(cls, tab_delimited_file, assembly, species, config_info):
        dt = time.strftime("%Y%m%d", time.gmtime())
        tab_delimited_file.write(dt)
        tab_delimited_file.write('\n')
        tab_delimited_file.write('HGVS\tConsequence\tSymbol\tAlleles\tAllele of genes\tSymbol (text)')
        tab_delimited_file.write('\n')

    @classmethod
    def _variant_value_for_file(cls, variant, data_key, transform=None):
        value = variant.get(data_key)
        if value is None:
            return None
        if transform is None:
            return value
        return transform(value)

    @classmethod
    def _add_variant_to_vcf_file(cls, vcf_file, variant):
        info_map = OrderedDict()
        info_map['hgvs_nomenclature'] = cls._variant_value_for_file(variant, 'hgvsNomenclature')

        variant['geneLevelConsequence'] = []
        variant['transcriptLevelConsequence'] = []
        variant['geneImpact'] = []
        variant['transcriptImpact'] = []
        variant['geneSymbols'] = []
        variant['transcriptGFF3IDs'] = []
        variant['geneIDs'] = []
        variant['transcriptIDs'] = []
        for geneConsequence in variant['geneConsequences']:
            if geneConsequence['consequence'] is not None:
                variant['geneLevelConsequence'].append(geneConsequence['consequence'])
            else:
                variant['geneLevelConsequence'].append('')
            if geneConsequence['impact'] is not None:
                variant['geneImpact'].append(geneConsequence['impact'])
            else:
                variant['geneImpact'].append('')
            if geneConsequence['gene'] is not None:
                variant['geneIDs'].append(geneConsequence['gene'])
                variant['geneSymbols'].append(geneConsequence['geneSymbol'])
            else:
                variant['geneSymbols'].append('')

        for transcriptConsequence in variant['transcriptConsequences']:
            if transcriptConsequence['consequence'] is not None:
                variant['transcriptLevelConsequence'].append(transcriptConsequence['consequence'].replace(",", "|"))
            else:
                variant['transcriptLevelConsequence'].append('')
            if transcriptConsequence['impact'] is not None:
                variant['transcriptImpact'].append(transcriptConsequence['impact'])
            else:
                variant['transcriptImpact'].append('')
            if transcriptConsequence['transcript'] is not None:
                variant['transcriptIDs'].append(transcriptConsequence['transcript'])
                if transcriptConsequence['transcriptGFF3ID']:
                    variant['transcriptGFF3IDs'].append(transcriptConsequence['transcriptGFF3ID'])
                else:
                    variant['transcriptGFF3IDs'].append('')
            else:
                variant['transcriptGFF3IDs'].append('')

        if cls._variant_value_for_file(variant, 'geneLevelConsequence') is not None:
            info_map['geneLevelConsequence'] = ','.join(cls._variant_value_for_file(variant, 'geneLevelConsequence'))
        else:
            info_map['geneLevelConsequence'] = cls._variant_value_for_file(variant, 'geneLevelConsequence')

        if cls._variant_value_for_file(variant, 'transcriptLevelConsequence') is not None:
            info_map['transcriptLevelConsequence'] = ','.join(cls._variant_value_for_file(variant, 'transcriptLevelConsequence'))
        else:
            info_map['transcriptLevelConsequence'] = cls._variant_value_for_file(variant, 'transcriptLevelConsequence')

        if cls._variant_value_for_file(variant, 'geneLevelConsequence') is not None:
            info_map['geneImpact'] = ','.join(cls._variant_value_for_file(variant, 'geneImpact'))
        else:
            info_map['geneImpact'] = cls._variant_value_for_file(variant, 'geneImpact')

        if cls._variant_value_for_file(variant, 'geneLevelConsequence') is not None:
            info_map['transcriptImpact'] = ','.join(cls._variant_value_for_file(variant, 'transcriptImpact'))
        else:
            info_map['transcriptImpact'] = cls._variant_value_for_file(variant, 'transcriptImpact')

        for allele in variant['alleles']:
            if 'id' in variant:
                variant['alleles_ids'].append(allele['id'])
            else:
                variant['alleles'] = [allele['id']]
            if 'alleleSymbols' in variant:
                variant['alleleSymbols'].append(allele['symbol'])
            else:
                variant['alleleSymbols'] = [allele['symbol']]
            if 'alleleSymbolText' in variant:
                variant['alleleSymbolText'].append(allele['symbolText'])
            else:
                variant['alleleSymbolText'] = [allele['symbolText']]

        info_map['alleles'] = cls._variant_value_for_file(variant, 'allele_ids', transform=','.join)
        info_map['allele_symbols'] = cls._variant_value_for_file(variant, 'alleleSymbols', transform=','.join)
        info_map['allele_symbols_text'] = cls._variant_value_for_file(variant, 'alleleSymbolText', transform=','.join)
        info_map['soTerm'] = cls._variant_value_for_file(variant, 'soTerm')
        info_map['globalId'] = variant['globalId']

        if variant['geneIDs']:
            info_map['allele_of_gene_ids'] = cls._variant_value_for_file(variant, 'geneIDs', transform=','.join)
            info_map['allele_of_gene_symbols'] = cls._variant_value_for_file(variant, 'geneSymbols', transform=','.join)

        if variant['transcriptIDs']:
            info_map['allele_of_transcript_ids'] = cls._variant_value_for_file(variant, 'transcriptIDs', transform=','.join)
            info_map['allele_of_transcript_gff3_ids'] = cls._variant_value_for_file(variant, 'transcriptGFF3IDs', transform=','.join)

        if any(info_map.values()):
            info = ';'.join('{}="{}"'.format(k, v)
                            for (k, v) in info_map.items()
                            if v)
        else:
            info = cls.empty_value_marker
        vcf_file.write('\t'.join([variant['chromosome'],
                                  str(variant['POS']),
                                  info_map['hgvs_nomenclature'],
                                  variant['genomicReferenceSequence'],
                                  variant['genomicVariantSequence'],
                                  '.',
                                  '.',
                                  info]))
        vcf_file.write('\n')

    @classmethod
    def _add_variant_to_tab_file(cls, tab_delimited_file, variant):
        info_map = OrderedDict()
        info_map['hgvs_nomenclature'] = cls._variant_value_for_file(variant, 'hgvsNomenclature')
        if cls._variant_value_for_file(variant, 'geneLevelConsequence') is not None:
            info_map['geneLevelConsequence'] = ' '.join(cls._variant_value_for_file(variant, 'geneLevelConsequence').split('_'))
        else:
            info_map['geneLevelConsequence'] = cls._variant_value_for_file(variant, 'geneLevelConsequence')
        if cls._variant_value_for_file(variant, 'geneLevelConsequence') is not None:
            info_map['impact'] = ' '.join(cls._variant_value_for_file(variant, 'impact'))
        else:
            info_map['impact'] = cls._variant_value_for_file(variant, 'impact')

        info_map['alleleSymbol'] = cls._variant_value_for_file(variant, 'alleleSymbol')
        info_map['globalId'] = variant['globalId']
        info_map['alleles'] = cls._variant_value_for_file(variant, 'alleles', transform=', '.join)
        info_map['allele_of_genes'] = cls._variant_value_for_file(variant, 'alleleOfGenes', transform=', '.join)
        info_map['symbol_text'] = cls._variant_value_for_file(variant, 'symbolText')
        if any(info_map.values()):
            info = '\t'.join(v for (k, v) in info_map.items() if v)
        else:
            info = cls.empty_value_marker
        tab_delimited_file.write(info)
        tab_delimited_file.write('\n')

    def _consume_data_source(self):
        assembly_chr_variants = defaultdict(lambda: defaultdict(list))
        assembly_species = {}
        for variant in self.variants:
            assembly = variant['assembly'].replace('.', '').replace('_', '')
            chromosome = variant['chromosome']
            assembly_chr_variants[assembly][chromosome].append(variant)
            assembly_species[assembly] = variant['species']
        return (assembly_chr_variants, assembly_species)

    def _find_replace(self, string, iupac_codes):
        # is the item in the dict?
        for item in string:
            # iterate by keys
            if item in iupac_codes:
                # look up and replace
                string = string.replace(item, "<" + item + ">")
                # return updated string
        return string

    def _adjust_variant(self, variant):
        so_term = variant['soTerm']
        if variant['start'] is None:
            return None

        # from https://www.bioinformatics.org/sms/iupac.html
        iupac_to_vcf_ref_codes = {"R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}

        variant['genomicVariantSequence'] = self._find_replace(variant['genomicVariantSequence'],
                                                               iupac_to_vcf_ref_codes)

        if so_term == 'deletion':
            variant['POS'] = variant['start'] - 1
            if variant['genomicReferenceSequence'] == '':
                logger.warn('No reference sequence for variant Id: %r', variant['ID'])
                return None
            if variant['genomicVariantSequence'] == '':
                self._add_padded_base_to_variant(variant, 'deletion')
        elif so_term == 'insertion':
            if variant['genomicReferenceSequence'] != '':
                logger.warn('Insertion Variant reference sequence is populated'
                            'when it should not be in '
                            'variant ID: %r',
                            variant['globalId'])
                return None
            if variant['genomicVariantSequence'] == '':
                variant['genomicVariantSequence'] = '.'
            variant['POS'] = variant['start']
            self._add_padded_base_to_variant(variant, 'insertion')
        elif so_term in ['point_mutation', 'MNV']:
            variant['POS'] = variant['start']
        elif so_term == 'delins':
            if variant['genomicVariantSequence'] == '':
                variant['genomicVariantSequence'] = '.'
                variant['POS'] = variant['start'] - 1
                self._add_padded_base_to_variant(variant, 'delins')
            elif len(variant['genomicVariantSequence']) == len(variant['genomicReferenceSequence']):
                variant['POS'] = variant['start']
            else:
                variant['POS'] = variant['start'] - 1
                self._add_padded_base_to_variant(variant, 'delins')
        else:
            logger.fatal('New SoTerm that We need to add logic for: %r', so_term)
            return None
        return variant

    def generate_files(self, skip_chromosomes=(), upload_flag=False, tab_flag=False):
        (assembly_chr_variants, assembly_species) = self._consume_data_source()
        for (assembly, chromo_variants) in assembly_chr_variants.items():
            filename = assembly + '-' + self.config_info.config['RELEASE_VERSION'] + '.vcf'
            tab_filename = assembly + '-' + self.config_info.config['RELEASE_VERSION'] + '.txt'
            filepath = os.path.join(self.generated_files_folder, filename)
            filepath_tab = os.path.join(self.generated_files_folder, tab_filename)
            if not tab_flag:
                logger.info('Generating VCF File for assembly %r', assembly)
                with open(filepath, 'w') as vcf_file:
                    self._write_vcf_header(vcf_file, assembly,
                                           assembly_species[assembly],
                                           self.config_info)
                    for (chromosome, variants) in sorted(chromo_variants.items(), key=itemgetter(0)):
                        # print(variants)
                        if chromosome in skip_chromosomes:
                            logger.info('Skipping VCF file generation for chromosome %r', chromosome)
                            continue
                        adjust_varient = partial(self._adjust_variant)
                        adjusted_variants = filter(None, map(adjust_varient, variants))
                        for variant in sorted(adjusted_variants, key=itemgetter('POS')):
                            self._add_variant_to_vcf_file(vcf_file, variant)
            else:
                logger.info('Generating TAB delimited File for assembly %r', assembly)
                with open(filepath_tab, 'w') as tab_delimited_file:
                    self._write_tab_delimited_header(tab_delimited_file, assembly,
                                                     assembly_species[assembly],
                                                     self.config_info)
                    for (chromosome, variants) in sorted(chromo_variants.items(), key=itemgetter(0)):
                        if chromosome in skip_chromosomes:
                            logger.info('Skipping tab delimited file generation for chromosome %r', chromosome)
                            continue
                        adjust_varient = partial(self._adjust_variant)
                        adjusted_variants = filter(None, map(adjust_varient, variants))
                        for variant in sorted(adjusted_variants, key=itemgetter('POS')):
                            self._add_variant_to_tab_file(tab_delimited_file, variant)

            stdout, stderr, return_code = compress('bgzip -c ' + filepath + ' > ' + filepath + '.gz')
            if return_code == 0:
                logger.info(filepath + ' compressed successfully')
            else:
                logger.error(filepath + ' could not be compressed, please check')
            if upload_flag:
                logger.info("Submitting to FMS")
                process_name = "1"
                upload.upload_process(process_name, filename, self.generated_files_folder, 'VCF', assembly, self.config_info)
                upload.upload_process(process_name, filename + ".gz", self.generated_files_folder, 'VCF-GZ', assembly, self.config_info)
