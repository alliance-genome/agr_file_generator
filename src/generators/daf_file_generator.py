import os
import logging
from datetime import datetime
from time import gmtime, strftime

import json
import csv
import upload

logger = logging.getLogger(name=__name__)


class DafFileGenerator:

    file_header_template = """#########################################################################
#
# Disease Association
# Source: Alliance of Genome Resources (Alliance)
# Orthology Filter: Stringent
# TaxonIDs: {taxonIDs}
# Datebase Version: {databaseVersion}
# Date: {datetimeNow}
#
#########################################################################
"""

    def __init__(self, disease_associations, generated_files_folder, config_info, taxon_id_fms_subtype_map):
        self.disease_associations = disease_associations
        self.config_info = config_info
        self.taxon_id_fms_subtype_map = taxon_id_fms_subtype_map
        self.generated_files_folder = generated_files_folder

    @classmethod
    def _generate_header(cls, config_info, taxon_ids):
        return cls.file_header_template.format(taxonIDs=",".join(taxon_ids),
                                               datetimeNow=strftime("%Y-%m-%d %H:%M:%S", gmtime()),
                                               databaseVersion=config_info.config['RELEASE_VERSION'])

    def generate_file(self, upload_flag=False):
        fields = ["Taxon",
                  "SpeciesName",
                  "DBobjectType",
                  "DBObjectID",
                  "DBObjectSymbol",
                  #"InferredGeneAssociation",
                  #"GeneProductFormID",
                  #"AdditionalGeneticComponent",
                  #"ExperimentalConditions",
                  "AssociationType",
                  #"Qualifier",
                  "DOID",
                  "DOtermName",
                  "WithOrthologs",
                  "InferredFromID",
                  "InferredFromSymbol",
                  #"Modifier-AssociationType",
                  #"Modifier-Qualifier",
                  #"Modifier-Genetic",
                  #"Modifier-ExperimentalConditions",
                  "EvidenceCode",
                  "EvidenceCodeName",
                  #"genetic-sex",
                  "Reference",
                  "Date",
                  "Source"]

        processed_disease_associations = {}
        processed_disease_associations_tsv = {}
        for disease_association in self.disease_associations:
            for evidence in disease_association["evidence"]:
                if disease_association["objectType"][0] == "Feature":
                    db_object_type = "allele"
                elif disease_association["objectType"][0] == "AffectedGenomicModel":
                    db_object_type = "affected_genomic_model"
                else:
                    db_object_type = disease_association["objectType"][0].lower()


                pub_id = evidence["pubMedID"] if evidence["pubMedID"] else evidence["pubModID"]
                if pub_id is None:
                    pub_id = ""

                do_name = disease_association["DOtermName"] if disease_association["DOtermName"] else ""
                #inferred_gene_association = ""
                #if db_object_type == "gene":
                #    inferred_gene_association = disease_association["dbObjectID"]
                #elif db_object_type == "allele":
                #    inferred_gene_association = ",".join(disease_association["inferredGeneAssociation"])

                if evidence["evidenceCode"] is not None:
                    evidence_code = evidence["evidenceCode"]
                else:
                    evidence_code = ""

                if evidence["evidenceCodeName"] is not None:
                    evidence_code_name = evidence['evidenceCodeName']
                else:
                    evidence_code_name = ""

                # gene_product_form_id = ""
                # additional_genetic_component = ""
                # experimental_conditions = ""
                # qualifier = ""
                # modifier_association_type = ""
                # modifier_qualifier = ""
                # modifier_genetic = ""
                # modifier_experimental_conditions = ""
                # genetic_sex = ""

                if disease_association["associationType"] in ["implicated_via_orthology", "biomarker_via_orthology"] and len(disease_association["withOrthologs"]) == 0:
                    print(disease_association)
                    exit()
                    continue

                if disease_association["dateAssigned"] is None and disease_association["associationType"] in ["implicated_via_orthology",
                                                                                                              "biomarker_via_orthology"]:
                    date_str = strftime("%Y-%m-%d", gmtime())
                else:
                    date_str = disease_association["dateAssigned"]

                inferred_from_ids = []
                inferred_from_symbols = []
                for entity in disease_association["inferredFromEntities"]:
                    inferred_from_ids.append(entity["primaryKey"])
                    if "symbol" in entity:
                        inferred_from_symbols.append(entity["symbol"])
                    elif "name" in entity:
                        inferred_from_symbols.append(entity["name"])
                    else:
                        logger.info("infferred from node not handled" + entity["primaryKey"])

                taxon_id = disease_association["taxonId"]
                processed_association = dict(zip(fields, [taxon_id,
                                                          disease_association["speciesName"],
                                                          db_object_type,
                                                          disease_association["dbObjectID"],
                                                          disease_association["dbObjectSymbol"],
                                                          #inferred_gene_association,
                                                          #gene_product_form_id,
                                                          #additional_genetic_component,
                                                          #experimental_conditions,
                                                          disease_association["associationType"].lower(),
                                                          #qualifier,
                                                          disease_association["DOID"],
                                                          do_name,
                                                          disease_association["withOrthologs"],
                                                          inferred_from_ids,
                                                          inferred_from_symbols,
                                                          #modifier_association_type,
                                                          #modifier_qualifier,
                                                          #modifier_genetic,
                                                          #modifier_experimental_conditions,
                                                          evidence_code,
                                                          evidence_code_name,
                                                          #genetic_sex,
                                                          pub_id,
                                                          datetime.strptime(date_str, "%Y-%m-%d").strftime("%Y%m%d"),
                                                          disease_association["dataProvider"]]))
                processed_association_tsv = processed_association.copy()

                processed_association_tsv["WithOrthologs"] = "|".join(set(disease_association["withOrthologs"])) if len(disease_association["withOrthologs"]) > 0 else ""
                processed_association_tsv["InferredFromID"] = "|".join(set(processed_association_tsv["InferredFromID"])) if len(processed_association_tsv["InferredFromID"]) > 0 else ""
                processed_association_tsv["InferredFromSymbol"] = "|".join(set(processed_association_tsv["InferredFromSymbol"])) if len(processed_association_tsv["InferredFromSymbol"]) > 0 else ""

                if taxon_id in processed_disease_associations:
                    processed_disease_associations_tsv[taxon_id].append(processed_association_tsv)
                    processed_disease_associations[taxon_id].append(processed_association)
                else:
                    processed_disease_associations_tsv[taxon_id] = [processed_association_tsv]
                    processed_disease_associations[taxon_id] = [processed_association]

        file_basename = "agr-daf-" + self.config_info.config['RELEASE_VERSION']
        combined_file_basepath = os.path.join(self.generated_files_folder, file_basename + '.combined')

        combined_filepath_tsv = combined_file_basepath + '.tsv'
        combined_disease_file = open(combined_filepath_tsv, 'w')
        combined_disease_file.write(self._generate_header(self.config_info, set(processed_disease_associations_tsv.keys())))
        combined_tsv_writer = csv.DictWriter(combined_disease_file, delimiter='\t', fieldnames=fields, lineterminator="\n")
        combined_tsv_writer.writeheader()

        combined_filepath_json = combined_file_basepath + '.json'
        with open(combined_filepath_json, 'w') as outfile:
            json.dump(sum(processed_disease_associations.values(), []), outfile)

        for taxon_id in processed_disease_associations:
            combined_tsv_writer.writerows(processed_disease_associations_tsv[taxon_id])

            taxon_file_basepath = os.path.join(self.generated_files_folder, file_basename + '.' + taxon_id)
            taxon_filepath_json = taxon_file_basepath + '.json'
            with open(taxon_filepath_json, 'w') as f:
                json.dump(processed_disease_associations[taxon_id], f)

            taxon_filename_tsv = taxon_file_basepath + '.tsv'
            with open(taxon_filename_tsv, 'w') as f:
                f.write(self._generate_header(self.config_info, [taxon_id]))
                tsv_writer = csv.DictWriter(f, delimiter='\t', fieldnames=fields, lineterminator="\n")
                tsv_writer.writeheader()
                tsv_writer.writerows(processed_disease_associations_tsv[taxon_id])

        combined_disease_file.close()

        if upload_flag:
            logger.info("Submitting disease files to FMS")
            process_name = "1"
            upload.upload_process(process_name, combined_filepath_tsv, self.generated_files_folder, 'DISEASE-ALLIANCE', 'COMBINED', self.config_info)
            upload.upload_process(process_name, combined_filepath_json, self.generated_files_folder, 'DISEASE-ALLIANCE-JSON', 'COMBINED', self.config_info)
            for taxon_id in processed_disease_associations:
                 for file_extension in ['json', 'tsv']:
                     filename = file_basename + "." + taxon_id + '.' + file_extension
                     datatype = "DISEASE-ALLIANCE"
                     if file_extension == "json":
                          datatype += "-JSON"
                     upload.upload_process(process_name, filename, self.generated_files_folder, datatype, self.taxon_id_fms_subtype_map[taxon_id], self.config_info)
