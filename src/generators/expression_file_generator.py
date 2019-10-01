import os
import logging
import csv
from time import gmtime, strftime

import upload

logger = logging.getLogger(name=__name__)


class ExpressionFileGenerator:

    file_header_template = """#########################################################################
#
# Expression
# Source: Alliance of Genome Resources (Alliance)
# Datebase Version: {databaseVersion}
# Date: {datetimeNow}
#
#########################################################################
"""

    def __init__(self, expressions, generated_files_folder, config_info):
        self.expressions = expressions
        self.config_info = config_info
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _generate_header(cls, config_info):
        return cls.file_header_template.format(datetimeNow=strftime("%Y-%m-%d %H:%M:%S", gmtime()),
                                               databaseVersion=config_info.config['RELEASE_VERSION'])

    def generate_file(self, upload_flag=False):
        filename = 'agr-expression-' + self.config_info.config['RELEASE_VERSION'] + '.tsv'
        output_filepath = os.path.join(self.generated_files_folder, filename)
        expression_file = open(output_filepath,'w')
        expression_file.write(self._generate_header(self.config_info))

        columns = ['Species',
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

        tsv_writer = csv.DictWriter(expression_file, delimiter='\t', fieldnames=columns, lineterminator="\n")
        tsv_writer.writeheader()
        for expression in self.expressions:
            row = dict(zip(columns, [None] * len(columns)))
            row['Species'] = expression['species']['name']
            row['SpeciesID'] = expression['species']['primaryKey']
            row['GeneID'] = expression['gene']['primaryKey']
            row['GeneSymbol'] = expression['gene']['symbol']
            row['Location'] = expression['location']
            for term in expression['terms']:
                if 'CrossReference' in term.labels:
                    if row['Source']:
                        row['Source'] += ','
                        row['Source'] += term['crossRefCompleteUrl']  # according to spec should use globalCrossRefId
                    else:
                        row['Source'] = term['crossRefCompleteUrl']
                elif 'Publication' in term.labels:
                    publication = term['pubMedId'] or term['pubModId']
                    # reference = row['Reference']
                    if row['Reference']:
                        row['Reference'] += ','
                        row['Reference'] += publication
                    else:
                        row['Reference'] = publication
                elif 'Stage' in term.labels:
                    #row['StageID'] = term['primaryKey']
                    row['StageTerm'] = term['name']
                elif 'MMOTerm' in term.labels:
                    row['AssayID'] = term['primaryKey']
                    row['AssayTermName'] = term['name']
            for ontologyPath in expression['ontologyPaths']:
                if ontologyPath['edge'] == 'ANATOMICAL_STRUCTURE':
                    row['AnatomyTermID'] = ontologyPath['primaryKey']
                    row['AnatomyTermName'] = ontologyPath['name']
                elif ontologyPath['edge'] == 'CELLULAR_COMPONENT':
                    row['CellularComponentID'] = ontologyPath['primaryKey']
                    row['CellularComponentTerm'] = ontologyPath['name']
                elif ontologyPath['edge'] == 'ANATOMICAL_SUB_SUBSTRUCTURE':
                    row['SubStructureID'] = ontologyPath['primaryKey']
                    row['SubStructureName'] = ontologyPath['name']
                elif ontologyPath['edge'] == 'CELLULAR_COMPONENT_QUALIFIER':
                    if row['CellularComponentQualifierIDs']:
                        row['CellularComponentQualifierIDs'] += ','
                        row['CellularComponentQualifierIDs'] += ontologyPath['primaryKey']
                    else:
                        row['CellularComponentQualifierIDs'] = ontologyPath['primaryKey']
                    if row['CellularComponentQualifierTermNames']:
                        row['CellularComponentQualifierTermNames'] += ','
                        row['CellularComponentQualifierTermNames'] += ontologyPath['name']
                    else:
                        row['CellularComponentQualifierTermNames'] = ontologyPath['name']
                elif ontologyPath['edge'] == 'ANATOMICAL_SUB_STRUCTURE_QUALIFIER':
                    if row['SubStructureQualifierIDs']:
                        row['SubStructureQualifierIDs'] += ','
                        row['SubStructureQualifierIDs'] += ontologyPath['primaryKey']
                    else:
                        row['SubStructureQualifierIDs'] = ontologyPath['primaryKey']
                    if row['SubStructureQualifierTermNames']:
                        row['SubStructureQualifierTermNames'] += ','
                        row['SubStructureQualifierTermNames'] += ontologyPath['name']
                    else:
                        row['SubStructureQualifierTermNames'] = ontologyPath['name']
                elif ontologyPath['edge'] == 'ANATOMICAL_STRUCTURE_QUALIFIER':
                    if row['AnatomyTermQualifierIDs']:
                        row['AnatomyTermQualifierIDs'] += ','
                        row['AnatomyTermQualifierIDs'] += ontologyPath['primaryKey']
                    else:
                        row['AnatomyTermQualifierIDs'] = ontologyPath['primaryKey']
                    if row['AnatomyTermQualifierTermNames']:
                        row['AnatomyTermQualifierTermNames'] += ','
                        row['AnatomyTermQualifierTermNames'] += ontologyPath['name']
                    else:
                        row['AnatomyTermQualifierTermNames'] = ontologyPath['name']
            tsv_writer.writerows([row])
        expression_file.close()
        if upload_flag:
            logger.info("Submitting to FMS")
            process_name = "1"
            upload.upload_process(process_name, filename, self.generated_files_folder, 'EXPRESSION', 'ALL', self.config_info)
