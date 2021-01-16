import os
import time
import json
import sys
from collections import defaultdict, OrderedDict
from functools import partial
from operator import itemgetter
from common import run_command
from validators import vcf_validator
import logging
import upload

from .vcf_file_generator import VcfFileGenerator

sys.path.append('../')

logger = logging.getLogger(name=__name__)


class VariantAlleleFileGenerator:

    def __init__(self, variant_alleles, generated_files_folder, config_info):
        self.variant_alleles = variant_alleles
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder

    def generate_files(self, skip_chromosomes=(), upload_flag=False, validate_flag=False):
        filename = 'variant-allele-' + self.config_info.config['RELEASE_VERSION']
        filepath = os.path.join(self.generated_files_folder, filename)
        logger.info('Generating VARIANT_ALLELE File')

        documents = []
        for variant_allele in self.variant_alleles:
            vcf_generator = VcfFileGenerator(self.variant_alleles, self.generated_files_folder, self.config_info)
            VcfFileGenerator._adjust_variant(vcf_generator, variant_allele)
            allele_associated_gene_ids = []
            allele_associated_gene_symbols = []
            if variant_allele['alleleAssociatedGenes'] is not None:
                for allele_associated_genes in variant_allele['alleleAssociatedGenes']:
                    if allele_associated_genes['id']:
                        allele_associated_gene_ids.append(allele_associated_genes['id'])
                    if allele_associated_genes['symbol']:
                        allele_associated_gene_symbols.append(allele_associated_genes['symbol'])

            variant_affected_gene_ids = []
            variant_affected_gene_symbols = []
            if variant_allele['variantAffectedGenes']:
                for variant_affected_genes in variant_allele['variantAffectedGenes']:
                    if variant_affected_genes['id']:
                        variant_affected_gene_ids.append(variant_affected_genes['id'])
                    if variant_affected_genes['symbol']:
                        variant_affected_gene_symbols.append(variant_affected_genes['symbol'])

            variant_symbol = variant_allele['variantId']
            if variant_symbol:
                variant_id_parts = variant_allele['variantId'].split(':g.')
                if len(variant_id_parts) == 2:
                    variant_symbol = variant_allele['chromosome'] + ":" + variant_id_parts[1]

            has_disease = "-"
            if variant_allele['alleleDiseaseCount'] + variant_allele['variantDiseaseCount'] > 0:
                has_disease = "yes"

            has_phenotype = "-"
            if variant_allele['allelePhenotypeCount'] + variant_allele['variantPhenotypeCount'] > 0:
                has_phenotype = "yes"

            category = "allele"
            if variant_allele['alleleVariantCount'] > 0:
                category = category + " with " + str(variant_allele['alleleVariantCount']) + " known variants"

            document = {'species': variant_allele['species'],
                        'species_id': variant_allele['taxonId'],
                        'allele_id': variant_allele['allele']['id'] if variant_allele['allele'] else None,
                        'allele_symbol': variant_allele['allele']['symbol'] if variant_allele['allele'] else None,
                        'allele_synonyms': variant_allele['alleleSyns'],
                        'variant_id': variant_allele['variantId'],
                        'variant_symbol': variant_symbol,
                        'variant_synonyms': variant_allele['variantSyns'],
                        'variant_cross_reference': variant_allele['variantCrossReferences'],
                        'allele_associated_gene_id': allele_associated_gene_ids,
                        'allele_associated_gene_symbol': allele_associated_gene_symbols,
                        'variant_affected_gene_id': variant_affected_gene_ids,
                        'variant_affected_gene_symbol': variant_affected_gene_symbols,
                        'category': category,
                        'variants_type_id': variant_allele['variationType']['id'] if variant_allele['variationType'] else None,
                        'variants_type_name': variant_allele['variationType']['name'] if variant_allele['variationType'] else None,
                        'variants_hgvs_names': variant_allele['hgvsNomenclature'],
                        'assembly': variant_allele['assembly'],
                        'chromosome': variant_allele['chromosome'],
                        'start_postiion': variant_allele['start'],
                        'end_postiion': variant_allele['end'],
                        'sequence_of_reference': variant_allele['genomicReferenceSequence'],
                        'sequence_of_variant': variant_allele['genomicVariantSequence'],
                        'most_severe_consequence_name': variant_allele['geneConsequences'],
                        'variant_information_reference': variant_allele['pubIds'],
                        'has_disease_annotations': has_disease,
                        'has_phenotype_annotations': has_phenotype}
            documents.append(document)

        if validate_flag:
            process_name = "1"
            filepath = os.path.join(self.generated_files_folder, filename)
            if upload_flag:
                logger.info("Submitting to FMS")
                upload.upload_process(process_name, filename + '.tsv', self.generated_files_folder, 'VARIANT-ALLELE', 'TSV', self.config_info)
                upload.upload_process(process_name, filepath + ".json", self.generated_files_folder, 'VARIANT-ALLLELE-JSON', 'JSON', self.config_info)
