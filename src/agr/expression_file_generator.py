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
        output_filepath = self.generated_files_folder + "/agr-expression-" + self.database_version + ".tsv"
        expression_file = open(output_filepath,'w')
        expression_file.write(self._generate_header(self.database_version))
    
        columns = ["Species",
                   "SpeciesId",
                   "GeneID",
                   "GeneSymbol",
                   "Location",
                   "StageID",
                   "StageTerm",
                   "AssayID",
                   "AssayTermName",
                   "CellularComponentID",
                   "CellularComponentTerm",
                   "CellularComponentQualifierIDs",
                   "CellularComponentQualifierTermNames",
                   "CellTypeId",
                   "CellTypeName",
                   "CellTypeQualifierIDs",
                   "CellTypeQualifierTermNames",
                   "AnatomyTermID",
                   "AnatomyTermName",
                   "AnatomyTermQaulifierIDs",
                   "AnatomyTermQualifierTermNames",
                   "Source",
                   "Reference"]

        expression_file.write("\t".join(columns) + "\n")

        cellular_component_qualifier_ids = ""
        cellular_component_qualifier_term_names = ""
        cell_type_id = ""
        cell_type_name = ""
        cell_type_qualifier_ids = ""
        cell_type_qualifier_term_names = ""
        anatomy_term_qualifier_ids = ""
        anatomy_term_qualifier_term_names = ""

        for expression in self.expressions:
            ref_obj = expression["Reference"]
            reference = ref_obj["pubMedId"] or ref_obj["pubModId"]
            anatomy_term_id = expression["AnatomyTermID"] or ""
            anatomy_term_name = expression["AnatomyTermName"] or ""
            cellular_component_id = expression["CellularComponentID"] or ""
            cellular_component_term = expression["CellularComponentTerm"] or ""

            expression_file.write("\t".join([expression["Species"],
                                          expression["SpeciesID"],
                                          expression["GeneID"],
                                          expression["GeneSymbol"],
                                          expression["location"],
                                          expression["StageID"],
                                          expression["StageTerm"],
                                          expression["AssayID"],
                                          expression["AssayTerm"],
                                          cellular_component_id,
                                          cellular_component_term,
                                          cellular_component_qualifier_ids,
                                          cellular_component_qualifier_term_names,
                                          cell_type_id,
                                          cell_type_name,
                                          cell_type_qualifier_ids,
                                          cell_type_qualifier_term_names,
                                          anatomy_term_id,
                                          anatomy_term_name,
                                          anatomy_term_qualifier_ids,
                                          anatomy_term_qualifier_term_names,
                                          expression["Source"],
                                          reference])
                             + "\n")
    
        expression_file.close()
