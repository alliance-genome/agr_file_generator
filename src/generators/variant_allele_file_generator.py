import os
import sys
import json
import csv

import logging
import upload

from headers import create_header
from .vcf_file_generator import VcfFileGenerator
from validators import json_validator

sys.path.append('../')

logger = logging.getLogger(name=__name__)


class VariantAlleleFileGenerator:

    def __init__(self, variant_alleles, generated_files_folder, config_info):
        self.variant_alleles = variant_alleles
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _generate_header(cls, config_info, taxon_ids, data_format):
        """
          :param config_info:
          :param taxon_ids:
          :return:
        """

        return create_header('Variant/Allele',
                             config_info.config['RELEASE_VERSION'],
                             taxon_ids=taxon_ids,
                             config_info=config_info,
                             data_format=data_format)

    def generate_files(self, skip_chromosomes=(), upload_flag=False, validate_flag=False):

        fields = ['Taxon',
                  'SpeciesName',
                  'AlleleId',
                  'AlleleSymbol',
                  'AlleleSynonyms',
                  'VariantId',
                  'VariantSymbol',
                  'VariantSynonyms',
                  'VariantCrossReferences',
                  'AlleleAssociatedGeneId',
                  'AlleleAssociatedGeneSymbol',
                  'VariantAffectedGeneId',
                  'VariantAffectedGeneSymbol',
                  'Category',
                  'VariantsTypeId',
                  'VariantsTypeName',
                  'VariantsHgvsNames',
                  'Assembly',
                  'Chromosome',
                  'StartPosition',
                  'EndPosition',
                  'SequenceOfReference',
                  'SequenceOfVariant',
                  'MostSevereConsequenceName',
                  'VariantInformationReference',
                  'HasDiseaseAnnotations',
                  'HasPhenotypeAnnotations']

        filename = 'variant-allele-' + self.config_info.config['RELEASE_VERSION']
        filepath_stub = os.path.join(self.generated_files_folder, filename)
        filepath_json = filepath_stub + ".json"
        filepath_tsv = filepath_stub + ".tsv"

        logger.info('Generating VARIANT_ALLELE File')

        taxon_ids = set()
        documents = []
        for variant_allele in self.variant_alleles:
            vcf_generator = VcfFileGenerator(self.variant_alleles, self.generated_files_folder, self.config_info)
            if variant_allele['start']:
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

            taxon_ids.add(variant_allele['taxonId'])

            document = {'Taxon': variant_allele['taxonId'],
                        'SpeciesName': variant_allele['species'],
                        'AlleleId': variant_allele['allele']['id'] if variant_allele['allele'] else None,
                        'AlleleSymbol': variant_allele['allele']['symbol'] if variant_allele['allele'] else None,
                        'AlleleSynonyms': variant_allele['alleleSyns'],
                        'VariantId': variant_allele['variantId'],
                        'VariantSymbol': variant_symbol,
                        'VariantSynonyms': variant_allele['variantSyns'],
                        'VariantCrossReferences': variant_allele['variantCrossReferences'],
                        'AlleleAssociatedGeneId': allele_associated_gene_ids,
                        'AlleleAssociatedGeneSymbol': allele_associated_gene_symbols,
                        'VariantAffectedGeneId': variant_affected_gene_ids,
                        'VariantAffectedGeneSymbol': variant_affected_gene_symbols,
                        'Category': category,
                        'VariantsTypeId': variant_allele['variationType']['id'] if variant_allele['variationType'] else None,
                        'VariantsTypeName': variant_allele['variationType']['name'] if variant_allele['variationType'] else None,
                        'VariantsHgvsNames': variant_allele['hgvsNomenclature'],
                        'Assembly': variant_allele['assembly'],
                        'Chromosome': variant_allele['chromosome'],
                        'StartPosition': variant_allele['start'],
                        'EndPosition': variant_allele['end'],
                        'SequenceOfReference': variant_allele['genomicReferenceSequence'],
                        'SequenceOfVariant': variant_allele['genomicVariantSequence'],
                        'MostSevereConsequenceName': variant_allele['geneConsequences'],
                        'VariantInformationReference': variant_allele['pubIds'],
                        'HasDiseaseAnnotations': has_disease,
                        'HasPhenotypeAnnotations': has_phenotype}
            documents.append(document)

        with open(filepath_json, 'w') as f:
            contents = {'metadata': self._generate_header(self.config_info, taxon_ids, 'json'),
                        'data': documents}
            json.dump(contents, f)

        for document in documents:
            for key in document:
                if isinstance(document[key], list):
                    document[key] = ','.join(document[key])

        with open(filepath_tsv, 'w') as f:
            f.write(self._generate_header(self.config_info, taxon_ids, 'tsv'))
            writer = csv.DictWriter(f, delimiter='\t', fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(documents)

        if validate_flag:
            process_name = "1"
            logger.info("validating JSON file")
            json_validator.JsonValidator(filepath_json, 'variant-allele').validateJSON()
            if upload_flag:
                logger.info("Submitting to FMS")
                upload.upload_process(process_name, filename + ".tsv", self.generated_files_folder, 'VARIANT-ALLELE', 'COMBINED', self.config_info)
                upload.upload_process(process_name, filename + ".json", self.generated_files_folder, 'VARIANT-ALLLELE-JSON', 'COMBINED', self.config_info)
