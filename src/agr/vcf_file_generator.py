from collections import defaultdict, OrderedDict
from functools import partial
from operator import itemgetter
import logging
import os
import time

import pyfaidx

import agr.assembly_sequence as agr_asm_seq


logger = logging.getLogger(name=__name__)


class VcfFileGenerator:

    empty_value_marker = '.'

    file_header = """##fileformat=VCFv4.2
##fileDate={datetime}
##source=agr_file_genrator
##reference=
##contig=<ID=,length=,assembly={assembly},md5=,species="{species}",taxonomy=x>
##phasing=partial
##INFO=<ID=hgvs_nomenclature,Type=String,Number=0,,Description="the HGVS name of the allele">
##INFO=<ID=symbol,Type=String,Number=0,Description="The human readable name of the allele">
##INFO=<ID=allele_of_genes,Type=String,Number=0,Description="The genes that the Allele is located on">
##INFO=<ID=DP,Number=0,Type=Integer,Description="The label to be used for visual purposes">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">"""

    col_headers = ('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')

    def __init__(self, variants, generated_files_folder, fasta_sequences_folder, database_version):
        self.variants = variants
        self.database_version = database_version
        self.fasta_sequences_folder = fasta_sequences_folder
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _add_genomic_reference_sequence(cls, fasta, variant):
        g_ref_seq = agr_asm_seq.get_seq(fasta, variant['chromosome'], variant['start'], variant['end'])
        variant['genomicReferenceSequence'] = g_ref_seq

    @classmethod
    def _add_padded_base_to_variant(cls, fasta, variant, so_term):
        padded_base = agr_asm_seq.get_seq(fasta, variant['chromosome'], variant['POS'], variant['POS'])
        variant['genomicReferenceSequence'] = padded_base + variant['genomicReferenceSequence']
        variant['genomicVariantSequence'] = padded_base + variant['genomicVariantSequence']

    @classmethod
    def _write_vcf_header(cls, vcf_file, assembly, species, database_version):
        dt = time.strftime("%Y%m%d", time.gmtime())
        header = cls.file_header.format(datetime=dt,
                                        database_version=database_version,
                                        species=species,
                                        assembly=assembly)
        vcf_file.write(header)
        vcf_file.write('\n')
        vcf_file.write('\t'.join(cls.col_headers))
        vcf_file.write('\n')

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
        info_map['symbol'] = cls._variant_value_for_file(variant, 'symbol')
        info_map['allele_of_genes'] = cls._variant_value_for_file(variant,
                                                                  'alleleOfGenes',
                                                                  transform=', '.join)
        if any(info_map.values()):
            info = ';'.join('{}="{}"'.format(k, v)
                            for (k, v) in info_map.items()
                            if v)
        else:
            info = cls.empty_value_marker
        vcf_file.write('\t'.join([variant['chromosome'],
                                  str(variant['POS']),
                                  variant['globalId'],
                                  variant['genomicReferenceSequence'],
                                  variant['genomicVariantSequence'],
                                  '.',
                                  '.',
                                  info]))
        vcf_file.write('\n')

    def _consume_data_source(self):
        assembly_chr_variants = defaultdict(lambda : defaultdict(list))
        assembly_species = {}
        for variant in self.variants:
            assembly = variant['assembly']
            chromosome = variant['chromosome']
            assembly_chr_variants[assembly][chromosome].append(variant)
            assembly_species[assembly] = variant['species']
        return (assembly_chr_variants, assembly_species)

    def _adjust_variant(self, fasta, variant):
        so_term = variant['soTerm']
        start_pos = variant['start']
        variant['POS'] = start_pos - 1 if so_term == 'insertion' else start_pos
        if so_term == 'deletion':
            if variant['genomicReferenceSequence'] == '':
                self._add_genomic_reference_sequence(fasta, variant)
            if variant['genomicVariantSequence'] == '':
                self._add_padded_base_to_variant(fasta, variant, 'deletion')
        elif so_term == 'insertion':
            if variant['genomicReferenceSequence'] != '':
                logger.error('Insertion Variant reference sequence is populated '
                             'when it should not be in '
                             'variant ID: %r',
                             variant['ID'])
                return None
            if variant['genomicVariantSequence'] == '':
                return None
            self._add_padded_base_to_variant(fasta, variant, 'insertion')
        elif so_term == 'point_mutation':
            variant['POS'] = variant['start']
        elif so_term == 'MNV':
            variant['POS'] = variant['end']
        else:
            logger.fatal('New SoTerm that We need to add logic for: %r', so_term)
            return None
        return variant

    def generate_files(self, skip_chromosomes=()):
        (assembly_chr_variants, assembly_species) = self._consume_data_source()

        for (assembly, chromo_variants) in assembly_chr_variants.items():
            logger.info('Generating VCF File for assembly %r', assembly)
            filename = assembly + '-' + self.database_version  + '.vcf'
            filepath = os.path.join(self.generated_files_folder, filename)
            fasta = pyfaidx.Fasta(agr_asm_seq.get_local_path(self.fasta_sequences_folder,
                                                             assembly))
            with open(filepath, 'w') as vcf_file:
                self._write_vcf_header(vcf_file,
                                       assembly,
                                       assembly_species[assembly],
                                       self.database_version)
                for (chromosome, variants) in sorted(chromo_variants.items(),
                                                     key=itemgetter(0)):
                    if chromosome in skip_chromosomes:
                        logger.info('Skipping VCF file generation for chromosome %r',
                                    chromosome)
                        continue
                    adjust_varient = partial(self._adjust_variant, fasta)
                    adjusted_variants = filter(None, map(adjust_varient, variants))
                    for variant in sorted(adjusted_variants, key=itemgetter('POS')):
                        self._add_variant_to_vcf_file(vcf_file, variant)
