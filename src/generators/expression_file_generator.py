"""
.. module:: expression_file_generators
    :platform: any
    :synopsis: Module that generates the Expression files
.. moduleauthor:: AGR consrotium

"""

import os
import logging
import json
import csv
import upload
from .header import create_header

logger = logging.getLogger(name=__name__)


class ExpressionFileGenerator:
    """
    TBA
    """

    def __init__(self, expressions, generated_files_folder, config_info, taxon_id_fms_subtype_map):
        """

        :param expressions:
        :param generated_files_folder:
        :param config_info:
        :param taxon_id_fms_subtype_map:
        """
        self.expressions = expressions
        self.config_info = config_info
        self.taxon_id_fms_subtype_map = taxon_id_fms_subtype_map
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _generate_header(cls, config_info, species):
        """

        :param config_info:
        :param taxon_ids:
        :return:
        """

        if len(species.keys()) == 1:
            species_names = ''.join(list(species.values()))
            taxon_ids = '# TaxonIDs:' + ''.join(species.keys())
        else:
            taxon_ids = '# TaxonIDs: NCBITaxon:9606, NCBITaxon:10116, NCBITaxon:10090, NCBITaxon:7955, NCBITaxon:7227, NCBITaxon:6239, NCBITaxon:559292'
            species_names = 'Homo sapiens, Rattus norvegicus, Mus musculus, Danio rerio, Drosophila melanogaster, Caenorhabditis elegans, Saccharomyces cerevisiae'

        return create_header('Expression', config_info.config['RELEASE_VERSION'],
                             taxon_ids=taxon_ids,
                             species=species_names,
                             data_format='tsv')

    # 'StageID', currently don't have stage IDs in the database
    def generate_file(self, upload_flag=False):
        """

        :param upload_flag:
        :return:
        """
        fields = ['Species',
                  'SpeciesID',
                  'GeneID',
                  'GeneSymbol',
                  'Location',
                  'StageTerm',
                  'AssayID',
                  'AssayTermName',
                  'CellularComponentID',
                  'CellularComponentTerm',
                  'CellularComponentQualifierIDs',
                  'CellularComponentQualifierTermNames',
                  'SubStructureID',
                  'SubStructureName',
                  'SubStructureQualifierIDs',
                  'SubStructureQualifierTermNames',
                  'AnatomyTermID',
                  'AnatomyTermName',
                  'AnatomyTermQualifierIDs',
                  'AnatomyTermQualifierTermNames',
                  'SourceURL',
                  'Source',
                  'Reference']

        associations = {}
        species = {}
        for expression in self.expressions:
            association = dict(zip(fields, [None] * len(fields)))
            association['Species'] = expression['species']['name']
            association['Source'] = expression['gene']['dataProvider']
            association['SpeciesID'] = expression['species']['primaryKey']
            association['SpeciesID'] = association['SpeciesID']
            association['GeneID'] = expression['gene']['primaryKey']
            association['GeneSymbol'] = expression['gene']['symbol']
            association['Location'] = expression['location']
            for term in expression['terms']:
                if 'CrossReference' in term.labels:
                    if association['SourceURL']:
                        association['SourceURL'].append(term['crossRefCompleteUrl'])  # according to spec should use globalCrossRefId
                    else:
                        association['SourceURL'] = [term['crossRefCompleteUrl']]
                elif 'Publication' in term.labels:
                    publication = term['pubMedId'] or term['pubModId']
                    # reference = association['Reference']
                    if association['Reference']:
                        association['Reference'].append(publication)
                    else:
                        association['Reference'] = [publication]
                elif 'Stage' in term.labels:
                    # association['StageID'] = term['primaryKey']
                    association['StageTerm'] = term['name']
                elif 'MMOTerm' in term.labels:
                    association['AssayID'] = term['primaryKey']
                    association['AssayTermName'] = term['name']
            for ontology_path in expression['ontologyPaths']:
                if ontology_path['edge'] == 'ANATOMICAL_STRUCTURE':
                    association['AnatomyTermID'] = ontology_path['primaryKey']
                    association['AnatomyTermName'] = ontology_path['name']
                elif ontology_path['edge'] == 'CELLULAR_COMPONENT':
                    association['CellularComponentID'] = ontology_path['primaryKey']
                    association['CellularComponentTerm'] = ontology_path['name']
                elif ontology_path['edge'] == 'ANATOMICAL_SUB_SUBSTRUCTURE':
                    association['SubStructureID'] = ontology_path['primaryKey']
                    association['SubStructureName'] = ontology_path['name']
                elif ontology_path['edge'] == 'CELLULAR_COMPONENT_QUALIFIER':
                    if association['CellularComponentQualifierIDs']:
                        association['CellularComponentQualifierIDs'].append(ontology_path['primaryKey'])
                    else:
                        association['CellularComponentQualifierIDs'] = [ontology_path['primaryKey']]
                    if association['CellularComponentQualifierTermNames']:
                        association['CellularComponentQualifierTermNames'].append(ontology_path['name'])
                    else:
                        association['CellularComponentQualifierTermNames'] = [ontology_path['name']]
                elif ontology_path['edge'] == 'ANATOMICAL_SUB_STRUCTURE_QUALIFIER':
                    if association['SubStructureQualifierIDs']:
                        association['SubStructureQualifierIDs'].append(ontology_path['primaryKey'])
                    else:
                        association['SubStructureQualifierIDs'] = [ontology_path['primaryKey']]
                    if association['SubStructureQualifierTermNames']:
                        association['SubStructureQualifierTermNames'].append(ontology_path['name'])
                    else:
                        association['SubStructureQualifierTermNames'] = [ontology_path['name']]
                elif ontology_path['edge'] == 'ANATOMICAL_STRUCTURE_QUALIFIER':
                    if association['AnatomyTermQualifierIDs']:
                        association['AnatomyTermQualifierIDs'].append(ontology_path['primaryKey'])
                    else:
                        association['AnatomyTermQualifierIDs'] = [ontology_path['primaryKey']]
                    if association['AnatomyTermQualifierTermNames']:
                        association['AnatomyTermQualifierTermNames'].append(ontology_path['name'])
                    else:
                        association['AnatomyTermQualifierTermNames'] = [ontology_path['name']]
            taxon_id = association['SpeciesID']
            species[taxon_id] = association['Species']
            if taxon_id in associations:
                associations[taxon_id].append(association)
            else:
                associations[taxon_id] = [association]

        file_basename = "agr-expression-" + self.config_info.config['RELEASE_VERSION']
        combined_file_basepath = os.path.join(self.generated_files_folder, file_basename + '.combined')

        combined_filepath_tsv = combined_file_basepath + '.tsv'
        combined_expression_file = open(combined_filepath_tsv, 'w')
        combined_expression_file.write(self._generate_header(self.config_info, species))
        combined_tsv_writer = csv.DictWriter(combined_expression_file, delimiter='\t', fieldnames=fields, lineterminator="\n")
        combined_tsv_writer.writeheader()

        combined_filepath_json = combined_file_basepath + '.json'
        with open(combined_filepath_json, 'w') as outfile:
            json.dump(sum(associations.values(), []), outfile)

        for taxon_id in associations:
            species_name = species[taxon_id]
            taxon_file_basepath = os.path.join(self.generated_files_folder, file_basename + '.' + taxon_id)
            taxon_filepath_json = taxon_file_basepath + '.json'
            with open(taxon_filepath_json, 'w') as f:
                json.dump(associations[taxon_id], f)

            logger.info(taxon_id)
            for association in associations[taxon_id]:
                for key in association:
                    if isinstance(association[key], list):
                        association[key] = ','.join(association[key])

            combined_tsv_writer.writerows(associations[taxon_id])
            taxon_filename_tsv = taxon_file_basepath + '.tsv'
            with open(taxon_filename_tsv, 'w') as f:
                f.write(self._generate_header(self.config_info, {taxon_id: species_name}))
                tsv_writer = csv.DictWriter(f, delimiter='\t', fieldnames=fields, lineterminator="\n")
                tsv_writer.writeheader()
                tsv_writer.writerows(associations[taxon_id])

        combined_expression_file.close()

        if upload_flag:
            logger.info("Submitting expression files to FMS")
            process_name = "1"
            upload.upload_process(process_name, combined_filepath_tsv, self.generated_files_folder, 'EXPRESSION-ALLIANCE', 'COMBINED', self.config_info)
            upload.upload_process(process_name, combined_filepath_json, self.generated_files_folder, 'EXPRESSION-ALLIANCE-JSON', 'COMBINED', self.config_info)
            for taxon_id in associations:
                for file_extension in ['json', 'tsv']:
                    filename = file_basename + "." + taxon_id + '.' + file_extension
                    datatype = "EXPRESSION-ALLIANCE"
                    if file_extension == "json":
                        datatype += "-JSON"
                    upload.upload_process(process_name,
                                          filename,
                                          self.generated_files_folder,
                                          datatype,
                                          self.taxon_id_fms_subtype_map[taxon_id],
                                          self.config_info)
