import time

from neo4j import GraphDatabase

from agr.assembly_sequence import AssemblySequence


class VcfFileGenerator:

    def __init__(self, uri, generated_files_folder, database_version):
        self.driver = GraphDatabase.driver(uri)
        self.database_version = database_version
        self.generated_files_folder = generated_files_folder
        self.get_variants_query = """MATCH (s:Species)-[:FROM_SPECIES]-(a:Allele)-[:VARIATION]-(v:Variant)-[l:LOCATED_ON]-(c:Chromosome)
MATCH (v:Variant)-[:VARIATION_TYPE]-(st:SOTerm)
RETURN c.primaryKey AS chromosome,
       v.globalId AS globalId,
       v.genomicReferenceSequence AS genomicReferenceSequence,
       v.genomicVariantSequence AS genomicVariantSequence,
       v.hgvs_nomenclature AS hgvsNomenclature,
       v.dataProvider AS dataProvider,
       a.symbol AS symbol,
       l.start AS start,
       l.end AS end,
       l.assembly AS assembly,
       s.name AS species,
       st.nameKey AS soTerm
       """

    def generate_files(self):
        with self.driver.session() as session:
            with session.begin_transaction() as tx:
                assembly_chr_variant_dict = {}
                assembly_species_dict = {}
                for record in tx.run(self.get_variants_query):
                    variant = record.data()
                    assembly = variant['assembly']
                    chromosome = variant['chromosome']
                    if assembly not in assembly_chr_variant_dict:
                        assembly_chr_variant_dict[assembly] = {chromosome: [variant]}
                    elif chromosome not in assembly_chr_variant_dict[assembly]:
                        assembly_chr_variant_dict[assembly][chromosome] = [variant]
                    else:
                        assembly_chr_variant_dict[assembly][chromosome].append(variant)
                    assembly_species_dict[assembly] = variant['species']

                for assembly in assembly_chr_variant_dict:
                    print(assembly)
                    filename = assembly + '-' + self.database_version  + '.vcf'
                    filepath = self.generated_files_folder + '/' + filename
                    assembly_sequence = AssemblySequence(assembly)
                    vcf_file = open(filepath, 'w')
                    self._write_vcf_header(vcf_file, assembly, assembly_species_dict[assembly], self.database_version)
                    for chromosome in assembly_chr_variant_dict[assembly]:
                         if chromosome == 'Unmapped_Scaffold_8_D1580_D1567':
                             continue
                         for variant in assembly_chr_variant_dict[assembly][chromosome]:
                             if variant['soTerm'] == 'deletion':
                                 if variant['genomicReferenceSequence'] == '':
                                      self._add_genomic_reference_sequence(assembly_sequence, variant)
                                 if variant['genomicVariantSequence'] == '':
                                      self._add_padded_base_to_variant(assembly_sequence, variant, 'deletion')
                                 self._add_variant_to_vcf_file(vcf_file, variant)
                             elif variant['soTerm'] == 'insertion':
                                 if variant['genomicReferenceSequence'] != '':
                                     print('ERROR: Insertion Variant reference sequence is populated when it should not be')
                                     exit()
                                 if variant['genomicVariantSequence'] == '':
                                     continue
                                 else:
                                     variant['POS'] = variant['start']
                                     self._add_padded_base_to_variant(assembly_sequence, variant, 'insertion')
                                     self._add_variant_to_vcf_file(vcf_file, variant)
                             elif variant['soTerm'] == 'point_mutation':
                                 variant['POS'] = variant['start']
                                 self._add_variant_to_vcf_file(vcf_file, variant)
                             elif variant['soTerm'] == 'MNV':
                                 variant['POS'] = variant['end']
                                 self._add_variant_to_vcf_file(vcf_file, variant)
                             else:
                                 print('New SoTerm that We need to add logic for')
                                 print(variant['soTerm'])
                                 exit(1)
                                 variant['POS'] = variant['start']
                                 self._add_variant_to_vcf_file(vcf_file, variant)
                    vcf_file.close()

    @classmethod
    def _add_genomic_reference_sequence(cls, assembly_sequence, variant):
        variant['genomicReferenceSequence'] = assembly_sequence.get(variant['chromosome'], variant['start'], variant['end'])

    @classmethod
    def _add_padded_base_to_variant(cls, assembly_sequence, variant, soTerm):
        if soTerm == 'insertion':
            variant['POS'] = variant['start']
        else:
            variant['POS'] = variant['start'] - 1

        padded_base = assembly_sequence.get(variant['chromosome'], variant['POS'], variant['POS'])
        variant['genomicReferenceSequence'] = padded_base + variant['genomicReferenceSequence']
        variant['genomicVariantSequence'] = padded_base + variant['genomicVariantSequence']

    @classmethod
    def _write_vcf_header(cls, vcf_file, assembly, species, database_version):
        dt = time.strftime("%Y%m%d", time.gmtime())
        header = """##fileformat=VCFv4.2
##fileDate={datetime}
##source=agr_file_generator/src/app.py
##reference=
##contig=<ID=,length=,assembly={assembly},md5=,species="{species}",taxonomy=x>
##phasing=partial
##INFO=<ID=hgvs_nomenclature,Type=String,Number=0,,Description="the HGVS name of the allele">
##INFO=<ID=Symbol,Type=String,Number=0,Description="The human readable name of the allele">
##INFO=<ID=DP,Number=0,Type=Integer,Description="The label to be used for visual purposes">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
""".format(datetime=dt,
           database_version=database_version,
           species=species,
           assembly=assembly)
        vcf_file.write(header)
        vcf_file.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']))

    @classmethod
    def _variant_value_for_file(cls, variant, data_key, default='.'):
        if variant[data_key]:
            return variant[data_key]
        return default

    @classmethod
    def _add_variant_to_vcf_file(cls, vcf_file, variant):
        hgvs = cls._variant_value_for_file(variant, 'hgvsNomenclature')
        if hgvs == '.':
            info = '.'
        else:
            symbol = cls._variant_value_for_file(variant, 'symbol')
            info = 'hgvs_nomenclature=\'' + hgvs + '\';' + 'Symbol=\'' + symbol + '\''
        vcf_file.write('\n' + '\t'.join([variant['chromosome'],
                                         str(variant['POS']),
                                         variant['globalId'],
                                         variant['genomicReferenceSequence'],
                                         variant['genomicVariantSequence'],
                                         '.',
                                         info,
                                         '.']))
