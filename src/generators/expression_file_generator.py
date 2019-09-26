
import logging
import csv
from time import gmtime, strftime

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

    def __init__(self, expressions, generated_files_folder, database_version):
        self.expressions = expressions
        self.database_version = database_version
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _generate_header(cls, database_version):
        return cls.file_header_template.format(datetimeNow=strftime("%Y-%m-%d %H:%M:%S", gmtime()),
                                               databaseVersion=database_version)

    def generate_file(self):
        output_filepath = self.generated_files_folder + '/agr-expression-' + self.database_version + '.tsv'
        expression_file = open(output_filepath,'w')
        expression_file.write(self._generate_header(self.database_version))

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
