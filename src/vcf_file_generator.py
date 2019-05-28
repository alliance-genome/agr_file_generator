import sys, os

from dateutil.parser import parse
from datetime import datetime
from time import gmtime, strftime
import requests
from neo4j import GraphDatabase

class VcfFileGenerator(object):
    def __init__(self, uri, generated_files_folder, apollo_sequence_endpoint, database_version):
        self.driver = GraphDatabase.driver(uri)
        self.database_version = database_version
        self.generated_files_folder = generated_files_folder
        self.apollo_sequence_endpoint = apollo_sequence_endpoint
        self.get_variants_query = """
MATCH (s:Species)-[:FROM_SPECIES]-(:Allele)-[:VARIATION]-(v:Variant)-[l:LOCATED_ON]-(c:Chromosome)
RETURN c.primaryKey AS chromosome,
       v.globalId AS globalId,
       v.genomicReferenceSequence AS genomicReferenceSequence,
       v.genomicVariantSequence AS genomicVariantSequence,
       v.dataProvider AS dataProvider,
       l.start AS start,
       l.end AS end,
       l.assembly AS assembly,
       s.name AS species
    """

    def generateFiles(self):
        with self.driver.session() as session:
            with session.begin_transaction() as tx:
                assembly_variant_map = {}
                assembly_species_map = {}
                for record in tx.run(self.get_variants_query):
                    variant = record.data()
                    assembly = variant["assembly"]
                    if assembly not in assembly_variant_map:
                        assembly_variant_map[assembly] = [variant]
                    else:
                        assembly_variant_map[assembly].append(variant)
                    assembly_species_map[assembly] = variant["species"]
        
                for assembly in assembly_variant_map:
                    filename = assembly + "-" + self.database_version  + ".vcf"
                    filepath = self.generated_files_folder + "/" + filename
                    vcf_file = open(filepath,'w')
                    VcfFileGenerator.__write_vcf_header(vcf_file, assembly, assembly_species_map[assembly], self.database_version)
        
                    for variant in assembly_variant_map[assembly]:
                         if variant["genomicVariantSequence"] == "N/A": # this is an deletion
                             variant["genomicVariantSequence"] = ""
                             if variant["genomicReferenceSequence"] == "":
                                  VcfFileGenerator.__add_genomic_reference_sequence_from_apollo(self.apollo_sequence_endpoint, variant)
                             VcfFileGenerator.__add_padded_base_to_variant(self.apollo_sequence_endpoint, variant)
                             VcfFileGenerator.__add_variant_to_vcf_file(vcf_file, variant)
                         elif variant["genomicReferenceSequence"] in ["N/A", ""]:
                             VcfFileGenerator.__add_genomic_reference_sequence_from_apollo(self.apollo_sequence_endpoint, variant)
                             VcfFileGenerator.__add_padded_base_to_variant(self.apollo_sequence_endpoint, variant)
                             VcfFileGenerator.__add_variant_to_vcf_file(vcf_file, variant)
                         else:
                             variant["POS"] = variant["start"]
                             VcfFileGenerator.__add_variant_to_vcf_file(vcf_file, variant)
                    vcf_file.close()

    def __add_genomic_reference_sequence_from_apollo(apollo_sequence_endpoint, variant):
        end = variant["end"] + 1
        url = VcfFileGenerator.__get_sequence_url(apollo_sequence_endpoint,
                                                  variant["species"],
                                                  variant["chromosome"],
                                                  variant["start"],
                                                  end)
        response = requests.get(url)
        if response.status_code != 200:
            print(variant)
            print("ERROR3: " + url)
            exit()
            return 4
    
        sequence = response.text
        if not VcfFileGenerator.__valid_ref_sequence(sequence):
            print("ERROR: Did not return valid sequence (" + url + ")")
            variant["genomicReferenceSequence"] = sequence

    def __add_padded_base_to_variant(apollo_sequence_endpoint, variant):
        if variant["genomicVariantSequence"] not in ["", "N/A"]:
            #insertion
            position_left = variant["start"]
            position_right = position_left + 1
            url = VcfFileGenerator.__get_sequence_url(apollo_sequence_endpoint, variant["species"], variant["chromosome"], position_left, position_right)
            response = requests.get(url)
            if response.status_code != 200:
                print(variant)
                print("ERROR: " + url)
                return 1
            else:
                paddedBase = response.text
                if paddedBase not in ["C","T","G","A", "N"]:
                    print("ERROR: Padded base contains non valid character(s)")
                    print(paddedBase)
                    return 2
                else:
                    variant["POS"] = position_left - 1
                    variant["genomicReferenceSequence"] = paddedBase + variant["genomicReferenceSequence"]
                    variant["genomicVariantSequence"] = paddedBase + variant["genomicVariantSequence"]
                    return 0
        else:
            # Deletion
            position_left = variant["start"] - 2
            position_right = position_left + 1
            url = VcfFileGenerator.__get_sequence_url(apollo_sequence_endpoint, variant["species"], variant["chromosome"], position_left, position_right)
            response = requests.get(url)
            if response.status_code != 200:
                print(variant)
                print("ERROR3: " + url)
                exit()
                return 1
            else:
                paddedBase = response.text
                if paddedBase not in ["C","T","G","A", "N"]:
                    print("ERROR: Padded base contains non valid character(s)")
                    print(paddedBase)
                    return 2
                else:
                    variant["POS"] = position_left
                    variant["genomicReferenceSequence"] = paddedBase + variant["genomicReferenceSequence"]
                    variant["genomicVariantSequence"] = paddedBase + variant["genomicVariantSequence"]
                    return 0

    def __clean_sequence(sequence):
        if sequence in ["N/A"]:
            return ""
        else:
            return sequence
    
    def __valid_ref_sequence(sequence):
        '''Returns true if sequence is valid'''
        if all(c in "CTGAN" for c in sequence):
            return True
        else:
            return False
    
    def __valid_alt_sequence(sequence):
        '''Returns true if sequence is valid'''
        if all(c in "CTGARYMKSWHBVDN" for c in sequence):
            return True
        else:
            return False

    def __get_sequence_url(apollo_sequence_endpoint, assembly, chromosome, left_position, right_position):
        return apollo_sequence_endpoint + assembly + "/" + chromosome + ":" + str(left_position) + ".." + str(right_position) 
    
    def __write_vcf_header(vcf_file, assembly, species, database_version):
        datetime = strftime("%Y%m%d", gmtime())
        header = """##fileformat=VCFv4.2
##fileDate={datetime}
##source=agr_file_generator/bin/generate_vcf_files.py
##reference=
##contig=<ID=,length=,assembly={assembly},md5=,species="{species}",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
    """.format(datetime = datetime,
               database_version=database_version,
               species=species,
               assembly=assembly)
        vcf_file.write(header)
        vcf_file.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]))

    def __add_variant_to_vcf_file(vcf_file, variant):
        vcf_file.write("\n" + "\t".join([variant["chromosome"],
                                         str(variant["POS"]),
                                         variant["globalId"],
                                         variant["genomicReferenceSequence"],
                                         variant["genomicVariantSequence"],
                                         ".",
                                         ".",
                                         "."]))
