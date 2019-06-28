import sys

import logging
from dateutil.parser import parse
from datetime import datetime
from time import gmtime, strftime


logger = logging.getLogger(name=__name__)


class DafFileGenerator:

    file_header_template = """#########################################################################
#
# Disease Association Format (DAF)
# Source: Alliance of Genome Resources (Alliance)
# Filter: stringent
# Datebase Version: {databaseVersion}
# Date: {datetimeNow}
#
#########################################################################
"""

    def __init__(self, disease_associations, generated_files_folder, database_version):
        self.disease_associations = disease_associations
        self.database_version = database_version
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _generate_header(cls, database_version):
        return cls.file_header_template.format(datetimeNow=strftime("%Y-%m-%d %H:%M:%S", gmtime()),
                databaseVersion=database_version)

    def generate_file(self):
        outputFilepath = self.generated_files_folder + "/agr-daf-" + self.database_version + ".tsv"
        disease_file = open(outputFilepath,'w')
        disease_file.write(self._generate_header(self.database_version))
    
        columns = ["Taxon",
                   "TaxonName",
                   "DBobjectType",
                   "DBObjectID",
                   "DBObjectSymbol",
                   "InferredGeneAssociation",
                   "AssociationType",
                   "DOID",
                   "DOname",
                   "withOrthologs",
                   "EvidenceCode",
                   "Reference",
                   "AssignedBy"]
        disease_file.write("\t".join(columns) + "\n")
        for disease_association in self.disease_associations:
            dbObjectType = "allele" if disease_association["objectType"][0] == "Feature" else disease_association["objectType"][0].lower()
            pubID = disease_association["pubMedID"] if disease_association["pubMedID"] else disease_association["pubModID"]
            if pubID is None:
                pubID = ""
 
            withOrthologs = "|".join(set(disease_association["withOrthologs"])) if disease_association["withOrthologs"] else ""
            DOname = disease_association["DOname"] if disease_association["DOname"] else ""
            inferredGeneAssociation = ""
            if dbObjectType == "gene":
                inferredGeneAssociation = disease_association["dbObjectID"]
            elif dbObjectType == "allele":
                inferredGeneAssociation = ",".join(disease_association["inferredGeneAssociation"])

            if disease_association["evidenceCode"] is not None:
                evidenceCode = disease_association["evidenceCode"]
            else:
                evidenceCode = ""

            disease_file.write("\t".join([disease_association["taxonId"],
                                          disease_association["speciesName"],
                                          dbObjectType,
                                          disease_association["dbObjectID"],
                                          disease_association["dbObjectSymbol"],
                                          inferredGeneAssociation,
                                          disease_association["associationType"].lower(),
                                          disease_association["DOID"],
                                          DOname,
                                          withOrthologs,
                                          evidenceCode,
                                          pubID,
                                          disease_association["dataProvider"]])
                             + "\n")
    
        #parse(disease_association["dateProduced"]).strftime("%Y%m%d"),
        disease_file.close()
