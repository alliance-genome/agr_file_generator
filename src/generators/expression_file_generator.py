import os
import logging
import json
import csv
from time import gmtime, strftime

import upload

logger = logging.getLogger(name=__name__)


class ExpressionFileGenerator:

    file_header_template = """#########################################################################
#
# Expression
# Source: Alliance of Genome Resources (Alliance)
# TaxonIDs: {taxonIDs}
# Datebase Version: {databaseVersion}
# Date: {datetimeNow}
#
#########################################################################
"""

    def __init__(self, expressions, generated_files_folder, config_info, taxon_id_fms_subtype_map):
        self.expressions = expressions
        self.config_info = config_info
        self.taxon_id_fms_subtype_map = taxon_id_fms_subtype_map
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _generate_header(cls, config_info, taxon_ids):
        return cls.file_header_template.format(taxonIDs=",".join(taxon_ids),
                                               datetimeNow=strftime("%Y-%m-%d %H:%M:%S", gmtime()),
                                               databaseVersion=config_info.config['RELEASE_VERSION'])

    def generate_file(self, upload_flag=False):
        fields = ['Species',
                  'SpeciesID',
                  'GeneID',
                  'GeneSymbol',
                  'Location',
                 # 'StageID', currently don't have stage IDs in the database
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
                  'Source',
                  'Reference']

        associations = dict()
        for expression in self.expressions:
            association = dict(zip(fields, [None] * len(fields)))
            association['Species'] = expression['species']['name']
            association['SpeciesID'] = expression['species']['primaryKey']
            association['SpeciesID'] = association['SpeciesID'].replace("NCBITaxonId", "NCBI:txid")
            association['GeneID'] = expression['gene']['primaryKey']
            association['GeneSymbol'] = expression['gene']['symbol']
            association['Location'] = expression['location']
            for term in expression['terms']:
                if 'CrossReference' in term.labels:
                    if association['Source']:
                        association['Source'].append(term['crossRefCompleteUrl']) # according to spec should use globalCrossRefId
                    else:
                        association['Source'] = [term['crossRefCompleteUrl']]
                elif 'Publication' in term.labels:
                    publication = term['pubMedId'] or term['pubModId']
                    # reference = association['Reference']
                    if association['Reference']:
                        association['Reference'].append(publication)
                    else:
                        association['Reference'] = [publication]
                elif 'Stage' in term.labels:
                    #association['StageID'] = term['primaryKey']
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
                        association['AnatomyTermQualifierTermNames'].append(ontology_path['name'])
                taxon_id = association['SpeciesID']
                if taxon_id in associations:
                    associations[taxon_id].append(association)
                else:
                    associations[taxon_id] = [association]

        file_basename = "agr-expression-" + self.config_info.config['RELEASE_VERSION']
        combined_file_basepath = os.path.join(self.generated_files_folder, file_basename + '.combined')

        combined_filepath_tsv = combined_file_basepath + '.tsv'
        combined_expression_file = open(combined_filepath_tsv, 'w')
        combined_expression_file.write(self._generate_header(self.config_info, set(associations.keys())))
        combined_tsv_writer = csv.DictWriter(combined_expression_file, delimiter='\t', fieldnames=fields, lineterminator="\n")
        combined_tsv_writer.writeheader()

        combined_filepath_json = combined_file_basepath + '.json'
        with open(combined_filepath_json, 'w') as outfile:
            json.dump(sum(associations.values(), []), outfile)

        for taxon_id in associations:
            taxon_file_basepath = os.path.join(self.generated_files_folder, file_basename + '.' + taxon_id)
            taxon_filepath_json = taxon_file_basepath + '.json'
            with open(taxon_filepath_json, 'w') as f:
                json.dump(associations[taxon_id], f)

        for taxon_id in associations:
            for association in associations[taxon_id]:
                for key in association:
                    if isinstance(association[key], list):
                        association[key] = ','.join(association[key])

            combined_tsv_writer.writerows(associations[taxon_id])
            taxon_filename_tsv = taxon_file_basepath + '.tsv'
            with open(taxon_filename_tsv, 'w') as f:
                f.write(self._generate_header(self.config_info, [taxon_id]))
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
                    datatype = "Expression-ALLIANCE"
                    if file_extension is "json":
                        datatype += "-JSON"
                    upload.upload_process(process_name, filename, self.generated_files_folder, datatype, self.taxon_id_fms_subtype_map[taxon_id], self.config_info)
