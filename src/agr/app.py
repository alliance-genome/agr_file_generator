import os

from agr.vcf_file_generator import VcfFileGenerator
from agr.data_source import DataSource


host = os.environ.get('NEO4J_HOST', 'localhost')

port = int(os.environ.get('NEO4J_PORT', 7687))

alliance_db_version = os.environ.get('ALLIANCE_DATABASE_VERSION', 'test')

uri = "bolt://" + host + ":" + str(port)

def main(generated_files_folder='generated_files'):
    os.makedirs(generated_files_folder, exist_ok=True)
    variants_query = """MATCH (s:Species)-[:FROM_SPECIES]-(a:Allele)-[:VARIATION]-(v:Variant)-[l:LOCATED_ON]-(c:Chromosome)
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
    data_source = DataSource(uri, variants_query)
    gvf = VcfFileGenerator(data_source, generated_files_folder, alliance_db_version)
    gvf.generate_files()


if __name__ == '__main__':
    main()
