import yaml
import requests

from neo4j import GraphDatabase

class DataSource:

    def __init__(self, uri, query):
        self.uri = uri
        self.driver = GraphDatabase.driver(self.uri)
        self.query = query

    def __repr__(self):
        s = '\n'.join(['<' + self.__class__.__qualname__ + '({uri},', '{query})'])
        return s.format(**dict((k, repr(v)) for (k, v) in vars(self).items()))

    def __iter__(self):
        with self.driver.session() as session:
            with session.begin_transaction() as tx:
                for record in tx.run(self.query):
                    yield record.data()

    def get_data(self):
        with self.driver.session() as session:
            with session.begin_transaction() as tx:
                return list(tx.run(self.query))

    def get_taxonomy_subtype_map():
        """Get a map linking taxonomy IDs (keys) to FMS subtype names (values).

        These are stored in the species.yaml file in the agr_schemas repository.
        """
        taxonomy_subtype_map = dict()

        url = 'https://raw.githubusercontent.com/alliance-genome/agr_schemas/master/ingest/species/species.yaml'

        species_stream = requests.get(url, stream=True)

        if species_stream.encoding is None:
            species_stream.encoding = 'utf-8'

        yaml_list = yaml.load(species_stream.content, Loader=yaml.SafeLoader)
        for item in yaml_list:
            taxonomy_subtype_map[item['taxonId']] = item['fmsSubtypeName']
            # taxonomy_subtype_map[item['taxonId']] = item['primaryDataProvider']['dataProviderShortName']

        return taxonomy_subtype_map
