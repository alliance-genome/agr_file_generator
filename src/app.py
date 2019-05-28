from vcf_file_generator import VcfFileGenerator
import os

if "NEO4J_HOST" in os.environ:
    host = os.environ['NEO4J_HOST']
else:
    host = "localhost"

if "NEO4J_PORT" in os.environ:
    port = int(os.environ['NEO4J_PORT'])
else:
    port = 7687

uri = "bolt://" + host + ":" + str(port)

if "ALLIANCE_DB_VERSION" in os.environ:
    alliance_db_version = os.environ['ALLIANCE_DATABASE_VERSION']
else:
    alliance_db_version = "test"


if __name__ == '__main__':
    apollo_sequence_endpoint = "https://agr-apollo.berkeleybop.io/apollo/sequence/"
    generated_files_folder = "generated_files"
    gvf = VcfFileGenerator(uri, generated_files_folder, apollo_sequence_endpoint, alliance_db_version)
    gvf.generateFiles()
