import sys

import logging
from dateutil.parser import parse
from datetime import datetime
from time import gmtime, strftime

from neo4j.v1 import GraphDatabase

uri = "staging.alliancegenome.org:7687"
driver = GraphDatabase.driver(uri)

outputFile = sys.argv[1]
databaseVersion = sys.argv[2]

logger = logging.getLogger(name=__name__)

class DafFileGenerator:

    disease_file = open(outputFile,'w')

    datetimeNow = strftime("%Y-%m-%d %H:%M:%S", gmtime())

    header = """#########################################################################
#
# Disease Association Format (DAF)
# Source: Alliance of Genome Resources (Alliance)
# Filter: stringent
# Datebase Version: {databaseVersion}
# Date: {datetimeNow}
#
#########################################################################
""".format(datetimeNow = datetimeNow, databaseVersion=databaseVersion)


    def createDaf:
        disease_file.write(header)
    
        disease_file.write("Taxon\tTaxonName\tDBobjectType\tDBObjectID\tDBObjectSymbol\tInferredGeneAssociation\tAssociationType\tDOID\tDOname\twithOrthologs\tEvidenceCode\tReference\tAssignedBy\n")
        with session.begin_transaction() as tx:
            for record in tx.run("""MATCH (dej:DiseaseEntityJoin)-[]-(object)-[da]->(disease:DOTerm)
                                    WHERE da.uuid = dej.primaryKey
                                    MATCH (object)-[FROM_SPECIES]->(species:Species)
                                    OPTIONAL MATCH (ec:EvidenceCode)-[EVIDENCE]-(dej)
                                    OPTIONAL MATCH (p:Publication)-[ev:EVIDENCE]-(dej)
                                    OPTIONAL MATCH (object)-[o:ORTHOLOGOUS]-(oGene:Gene)
                                    WHERE o.strictFilter AND (ec.primaryKey = "ISS" OR ec.primaryKey = "ISO")
                                    OPTIONAL MATCH (object)-[IS_ALLELE_OF]->(gene:Gene)
                                    RETURN  object.taxonId AS taxonId,
                                            species.name AS speciesName,
                                            collect(oGene.primaryKey) AS withOrthologs,
                                            labels(object) AS objectType,
                                            object.primaryKey AS dbObjectID,
                                            object.symbol AS dbObjectSymbol,
                                            p.pubMedId AS pubMedID,
                                            p.pubModId As pubModID,
                                            type(da) AS associationType,
                                            collect(DISTINCT gene.primaryKey) AS inferredGeneAssociation,
                                            disease.doId AS DOID,
                                            disease.name as DOname,
                                            ec.primaryKey AS evidenceCode,
                                            da.dataProvider AS dataProvider
                                            """):
    
                dbObjectType = "allele" if record["objectType"][0] == "Feature" else record["objectType"][0].lower()
                pubID = record["pubMedID"] if record["pubMedID"] else record["pubModID"]
     
                withOrthologs = "|".join(set(record["withOrthologs"])) if record["withOrthologs"] else ""
                DOname = record["DOname"] if record["DOname"] else ""
                inferredGeneAssociation = ""
                if dbObjectType == "gene":
                    inferredGeneAssociation = record["dbObjectID"]
                elif dbObjectType == "allele":
                    inferredGeneAssociation = ",".join(record["inferredGeneAssociation"])
    
                disease_file.write("\t".join([record["taxonId"],
                                              record["speciesName"],
                                              dbObjectType,
                                              record["dbObjectID"],
                                              record["dbObjectSymbol"],
                                              inferredGeneAssociation,
                                              record["associationType"].lower(),
                                              record["DOID"],
                                              DOname,
                                              withOrthologs,
                                              record["evidenceCode"],
                                              pubID,
    #                                         parse(record["dateProduced"]).strftime("%Y%m%d"),
                                              record["dataProvider"]])
                                 + "\n")
    
        disease_file.close()
