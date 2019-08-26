import sys

import logging
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

        expression_file.write('\t'.join(columns) + '\n')

        for expression in self.expressions:
            print(expression)
            anatomy_term_id = expression['anatomyTermObj']['id'] or ''
            anatomy_term_name = expression['anatomyTermObj']['term'] or ''
            anatomy_term_qualifier_ids = ','.join(str(anatomyTermQualifier['primaryKey']) for anatomyTermQualifier in expression['anatomyTermQualifiers']) or ''
            anatomy_term_qualifier_term_names = ','.join(str(anatomyTermQualifier['name']) for anatomyTermQualifier in expression['anatomyTermQualifiers']) or ''

            cellular_component_id = expression['cellularComponentObj']['id'] or ''
            cellular_component_term = expression['cellularComponentObj']['term'] or ''
            cellular_component_qualifier_ids = ','.join(str(cellularComponentQualifierObj['primaryKey']) for cellularComponentQualifierObj in expression['cellularComponentQualifiers']) or ''
            cellular_component_qualifier_term_names = ','.join(str(cellularComponentQualifierObj['name']) for cellularComponentQualifierObj in expression['cellularComponentQualifiers']) or ''

            cell_type_id = expression['cellTypeObj']['id'] or ''
            cell_type_term = expression['cellTypeObj']['term'] or ''
            cell_type_qualifier_ids = ','.join(str(cellTypeQualifierObj['primaryKey']) for cellTypeQualifierObj in expression['cellTypeQualifiers']) or ''
            cell_type_qualifier_term_names = ','.join(str(cellTypeQualifierObj['name']) for cellTypeQualifierObj in expression['cellTypeQualifiers']) or ''

            stage_id = expression['stageObj']['id'] or ''
            stage_term = expression['stageObj']['term'] or ''
            assay_id = expression['assayObj']['id'] or ''
            assay_term = expression['assayObj']['term'] or ''

            references = ','.join(ref_obj['pubMedId'] or ref_obj['pubModId'] for ref_obj in expression['References']) or ''

            expression_file.write('\t'.join([expression['speciesObj']['name'],
                                             expression['speciesObj']['id'],
                                             expression['geneObj']['id'],
                                             expression['geneObj']['symbol'],
                                             '\"' + expression['location'].replace('"', '\\"').replace(',', '\\,') + '\"',
                                             stage_id,
                                             stage_term,
                                             assay_id,
                                             assay_term,
                                             cellular_component_id,
                                             cellular_component_term,
                                             cellular_component_qualifier_ids,
                                             cellular_component_qualifier_term_names,
                                             cell_type_id,
                                             cell_type_term,
                                             cell_type_qualifier_ids,
                                             cell_type_qualifier_term_names,
                                             anatomy_term_id,
                                             anatomy_term_name,
                                             anatomy_term_qualifier_ids,
                                             anatomy_term_qualifier_term_names,
                                             expression['Source'],
                                             references])
                             + '\n')
    
        expression_file.close()
