from neo4j import GraphDatabase


class DataSource:

    def __init__(self, uri, query):
        self.uri = uri
        self.driver = GraphDatabase.driver(self.uri)
        self.query = query

    def __repr__(self):
        s =  '\n'.join(['<' + self.__class__.__qualname__ + '({uri},',
                        '{query})'])
        return s.format(**dict((k, repr(v)) for (k, v) in vars(self).items()))

    def __iter__(self):
        with self.driver.session() as session:
            with session.begin_transaction() as tx:
                for record in tx.run(self.query):
                    yield record.data()
