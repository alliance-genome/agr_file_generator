import os

from agr.vcf_file_generator import VcfFileGenerator


host = os.environ.get('NEO4J_HOST', 'localhost')

port = int(os.environ.get('NEO4J_PORT', 7687))

alliance_db_version = os.environ.get('ALLIANCE_DATABASE_VERSION', 'test')

uri = "bolt://" + host + ":" + str(port)

def main(generated_files_folder='generated_files'):
    os.makedirs(generated_files_folder, exist_ok=True)
    gvf = VcfFileGenerator(uri, generated_files_folder, alliance_db_version)
    gvf.generate_files()


if __name__ == '__main__':
    main()
