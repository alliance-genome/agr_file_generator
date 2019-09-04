import sys

import logging, csv
from dateutil.parser import parse
from datetime import datetime
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
                   'SpeciesId',
                   'GeneID',
                   'GeneSymbol',
                   'Location',
                   'StageID',
                   'StageTerm',
                   'AssayID',
                   'AssayTermName',
                   'CellularComponentID',
                   'CellularComponentTerm',
                   'CellularComponentQualifierIDs',
                   'CellularComponentQualifierTermNames',
                   'CellTypeId',
                   'CellTypeName',
                   'CellTypeQualifierIDs',
                   'CellTypeQualifierTermNames',
                   'AnatomyTermID',
                   'AnatomyTermName',
                   'AnatomyTermQaulifierIDs',
                   'AnatomyTermQualifierTermNames',
                   'Source',
                   'Reference']
        
        tsv_writer = csv.DictWriter(expression_file, delimiter='\t', fieldnames=columns, lineterminator="\n")
        tsv_writer.writeheader()
        for expression in self.expressions:
            publications = []
            sources = []
            stage_id = ''
            stage_term = ''
            assay_id = ''
            assay_term = ''
            anatomy_term_id = ''
            anatomy_term_name = ''
            for term in expression['terms']:
               if 'CrossReference' in term.labels:
                   sources.append(term['crossRefCompleteUrl']) # according to spec should use globalCrossRefId
               elif 'Publication' in term.labels:
                   publications.append(term['pubMedId'] or term['pubModId'])
               elif 'Stage' in term.labels:
                   stage_id = term['primaryKey']
                   stage_term = term['name']
               elif 'MMOTerm' in term.labels:
                   assay_id = term['primaryKey']
                   assay_term = term['name']
               elif 'UBERONTerm' in term.labels:
                   anatomy_term_id = term['primaryKey']
                   anatomy_term_name = term['name']
               else:
                   logger.error('Not handling term %r off gene %r %r', term['primaryKey'], expression['gene']['primaryKey'])

            location = ''
            cell_type_id = ''
            cell_type_name = ''
            cellular_component_id = ''
            cellular_component_term = ''
            cellular_type_qualifier_ids = []
            cellular_type_qualifier_term_names = []
            cellular_component_qualifier_ids = []
            cellular_component_qualifier_term_names = []
            anatomy_term_qualifier_ids = []
            anatomy_term_qualifier_term_names = []
            for [expressionBioEntity, association, ontology] in expression['entities']:
                location = expressionBioEntity['whereExpressedStatement']
                if association.type in ['CELLULAR_COMPONENT_RIBBON_TERM', 'ANATOMICAL_RIBBON_TERM']:
                    pass
                elif association.type == 'ANATOMICAL_STRUCTURE':
                    anatomy_term_id = ontology['primaryKey']
                    anatomy_term_name = ontology['name']
                elif association.type == 'CELLULAR_COMPONENT':
                    if 'GOTerm' in ontology.labels:
                        cellular_component_id = ontology['primaryKey']
                        cellular_component_term = ontology['name']
                    else:
                        print('cc entity')
                        print(association)
                        print(ontology)
                        exit()
                elif association.type == 'ANATOMICAL_SUB_SUBSTRUCTURE':
                    cell_type_id = ontology['primaryKey']
                    cell_type_name = ontology['name']
                elif association.type == 'CELLULAR_COMPONENT_QUALIFIER':
                    cellular_component_qualifier_ids.append(ontology['primaryKey'])
                    cellular_component_qualifier_term_names.append(ontology['name'])
                elif association.type == 'ANATOMICAL_SUB_STRUCTURE_QUALIFIER':
                    cell_type_qualifier_ids.append(ontology['primaryKey'])
                    cell_type_qualifier_term_names.append(ontology['name'])
                elif association.type == 'ANATOMICAL_STRUCTURE_QUALIFIER':
                    anatomy_term_qualifier_ids.append(ontology['primaryKey'])
                    anatomy_term_qualifier_term_names(ontology['primaryKey'])
                else:
                    logger.error('Not handling ontology term %r off gene %r %r', ontology['primaryKey'], expression['gene']['primaryKey'])

            row = dict(zip(columns, [expression['species']['primaryKey'],
                                     expression['species']['name'],
                                     expression['gene']['primaryKey'],
                                     expression['gene']['symbol'],
                                     location,
                                     stage_id,
                                     stage_term,
                                     assay_id,
                                     assay_term,
                                     cellular_component_id,
                                     cellular_component_term,
                                     ','.join(cellular_component_qualifier_ids) or '',
                                     ','.join(cellular_component_qualifier_term_names) or '',
                                     cell_type_id,
                                     cell_type_name,
                                     ','.join(cellular_type_qualifier_ids) or '',
                                     ','.join(cellular_type_qualifier_term_names) or '',
                                     anatomy_term_id or '',
                                     anatomy_term_name or '',
                                     ','.join(anatomy_term_qualifier_ids) or '',
                                     ','.join(anatomy_term_qualifier_term_names) or '',
                                     ','.join(list(dict.fromkeys(sources))) or '',
                                     ','.join(list(dict.fromkeys(publications))) or '']))
            tsv_writer.writerows([row])
        expression_file.close()
