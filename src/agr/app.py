import logging
import os

from agr.vcf_file_generator import VcfFileGenerator
import agr.assembly_sequence as agr_asm_seq
from agr.data_source import DataSource


host = os.environ.get('NEO4J_HOST', 'localhost')

port = int(os.environ.get('NEO4J_PORT', 7687))

alliance_db_version = os.environ.get('ALLIANCE_DATABASE_VERSION', 'test')

uri = "bolt://" + host + ":" + str(port)

assembly_to_s3_uri = {
    'R6.27': 'https://s3.amazonaws.com/agrjbrowse/FlyBase/fruitfly/fasta/',
    'GRCz11': 'https://s3.amazonaws.com/agrjbrowse/zfin/zebrafish/fasta/',
    'WBcel235': 'https://s3.amazonaws.com/agrjbrowse/WormBase/c_elegans_PRJNA13758/fasta/',
    'Rnor_6.0': 'https://s3.amazonaws.com/agrjbrowse/RGD/rat/fasta/',
    'GRCm38': 'https://s3.amazonaws.com/agrjbrowse/MGI/mouse/fasta/'
}


def setup_logging(logger_name):
    logging.basicConfig(level=logging.DEBUG)


def main(generated_files_folder='generated_files',
         fasta_sequences_folder='sequences',
         skip_chromosomes={'Unmapped_Scaffold_8_D1580_D1567'}):
    os.makedirs(generated_files_folder, exist_ok=True)
    os.makedirs(fasta_sequences_folder, exist_ok=True)
    variants_query = """MATCH (s:Species)-[:FROM_SPECIES]-(a:Allele)-[:VARIATION]-(v:Variant)-[l:LOCATED_ON]-(c:Chromosome)
MATCH (v:Variant)-[:VARIATION_TYPE]-(st:SOTerm)
OPTIONAL MATCH (a:Allele)-[:IS_ALLELE_OF]-(g:Gene)
RETURN c.primaryKey AS chromosome,
       v.globalId AS globalId,
       v.genomicReferenceSequence AS genomicReferenceSequence,
       v.genomicVariantSequence AS genomicVariantSequence,
       v.hgvs_nomenclature AS hgvsNomenclature,
       v.dataProvider AS dataProvider,
       a.symbol AS symbol,
       collect(g.primaryKey) as alleleOfGenes,
       l.start AS start,
       l.end AS end,
       l.assembly AS assembly,
       s.name AS species,
       st.nameKey AS soTerm
    """
    data_source = DataSource(uri, variants_query)
    agr_asm_seq.ensure_downloaded(fasta_sequences_folder, assembly_to_s3_uri)
    gvf = VcfFileGenerator(data_source,
                           generated_files_folder,
                           fasta_sequences_folder,
                           alliance_db_version)
    gvf.generate_files(skip_chromosomes=skip_chromosomes)


if __name__ == '__main__':
    main()
