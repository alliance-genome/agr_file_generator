import sys, os

from dateutil.parser import parse
from datetime import datetime
from time import gmtime, strftime
import requests
from neo4j import GraphDatabase
from assemblySequence import AssemblySequence

class VcfFileGenerator(object):
    def __init__(self, uri, generated_files_folder, database_version):
        self.driver = GraphDatabase.driver(uri)
        self.database_version = database_version
        self.generated_files_folder = generated_files_folder
        self.get_variants_query = """
MATCH (s:Species)-[:FROM_SPECIES]-(:Allele)-[:VARIATION]-(v:Variant)-[l:LOCATED_ON]-(c:Chromosome)
MATCH (v:Variant)-[:VARIATION_TYPE]-(st:SOTerm)
RETURN c.primaryKey AS chromosome,
       v.globalId AS globalId,
       v.genomicReferenceSequence AS genomicReferenceSequence,
       v.genomicVariantSequence AS genomicVariantSequence,
       v.dataProvider AS dataProvider,
       l.start AS start,
       l.end AS end,
       l.assembly AS assembly,
       s.name AS species,
       st.nameKey AS soTerm"""

    def generateFiles(self):
        with self.driver.session() as session:
            with session.begin_transaction() as tx:
                assembly_chr_variant_dict = {}
                assembly_species_dict = {}
                for record in tx.run(self.get_variants_query):
                    variant = record.data()
                    assembly = variant["assembly"]
                    chromosome = variant["chromosome"]
                    if assembly not in assembly_chr_variant_dict:
                        assembly_chr_variant_dict[assembly] = {chromosome: [variant]}
                    elif chromosome not in assembly_chr_variant_dict[assembly]:
                        assembly_chr_variant_dict[assembly][chromosome] = [variant]
                    else:
                        assembly_chr_variant_dict[assembly][chromosome].append(variant)
                    assembly_species_dict[assembly] = variant["species"]

                for assembly in assembly_chr_variant_dict:
                    filename = assembly + "-" + self.database_version  + ".vcf"
                    filepath = self.generated_files_folder + "/" + filename
                    assembly_sequence = AssemblySequence(assembly)
                    vcf_file = open(filepath,'w')
                    VcfFileGenerator.__write_vcf_header(vcf_file, assembly, assembly_species_dict[assembly], self.database_version)
                    for chromosome in assembly_chr_variant_dict[assembly]:
                         if chromosome == "Unmapped_Scaffold_8_D1580_D1567":
                             continue
                         for variant in assembly_chr_variant_dict[assembly][chromosome]:
                             if variant["soTerm"] == "deletion":
                                 if variant["genomicReferenceSequence"] == "":
                                      VcfFileGenerator.__add_genomic_reference_sequence(assembly_sequence, variant)
                                 if variant["genomicVariantSequence"] == "":
                                      VcfFileGenerator.__add_padded_base_to_variant(assembly_sequence, variant, "deletion")
                                 VcfFileGenerator.__add_variant_to_vcf_file(vcf_file, variant)
                             elif variant["soTerm"] == "insertion":
                                 if variant["genomicReferenceSequence"] != "":
                                     print("ERROR: Insertion Variant reference sequence is populated when it should not be")
                                     exit()
                                 if variant["genomicVariantSequence"] == "":
                                     continue
                                 else:
                                     variant["POS"] = variant["start"]
                                     VcfFileGenerator.__add_padded_base_to_variant(assembly_sequence, variant, "insertion")
                                     VcfFileGenerator.__add_variant_to_vcf_file(vcf_file, variant)
                             elif variant["soTerm"] == "point_mutation":
                                 variant["POS"] = variant["start"]
                                 VcfFileGenerator.__add_variant_to_vcf_file(vcf_file, variant)
                             elif variant["soTerm"] == "MNV":
                                 variant["POS"] = variant["end"]
                                 VcfFileGenerator.__add_variant_to_vcf_file(vcf_file, variant)
                             else:
                                 print("New SoTerm that We need to add logic for")
                                 print(variant["soTerm"])
                                 exit(1)
                                 variant["POS"] = variant["start"]
                                 VcfFileGenerator.__add_variant_to_vcf_file(vcf_file, variant)
                    vcf_file.close()

    def __add_genomic_reference_sequence(assembly_sequence, variant):
        variant["genomicReferenceSequence"] = assembly_sequence.get(variant["chromosome"], variant["start"], variant["end"])

    def __add_padded_base_to_variant(assembly_sequence, variant, soTerm):
        if soTerm == "insertion":
            variant["POS"] = variant["start"]
        else:
            variant["POS"] = variant["start"] - 1
        
        padded_base = assembly_sequence.get(variant["chromosome"], variant["POS"], variant["POS"])
        variant["genomicReferenceSequence"] = padded_base + variant["genomicReferenceSequence"] 
        variant["genomicVariantSequence"] = padded_base + variant["genomicVariantSequence"]

    def __write_vcf_header(vcf_file, assembly, species, database_version):
        datetime = strftime("%Y%m%d", gmtime())
        header = """##fileformat=VCFv4.2
##fileDate={datetime}
##source=agr_file_generator/src/app.py
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
