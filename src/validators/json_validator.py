import os
import logging
import json

from jsonschema import validate
from jsonschema import ValidationError
from jsonschema import SchemaError

logger = logging.getLogger(name=__name__)


class JsonValidator:
    def __init__(self, filepath, schema):
        self.filepath = filepath
        self.schema = schema

    def validateJSON(self):
        logger.info("validating " + self.filepath)
        with open(self.filepath, "r") as jsonFile:
            schema_filepath = os.path.join("./schemas/", self.schema) + '.schema'
            with open(schema_filepath, "r") as schemaFile:
                try:
                    validate(json.load(jsonFile), json.load(schemaFile))
                    logger.info("successfully validated against '%s'" % schema_filepath)
                except SchemaError as e:
                    logger.error(e)
                    logger.error("There is an error with the schema")
                    exit(-1)
                except ValidationError as e:
                    logger.error(e)
                    logger.error("---------")
                    logger.error(e.absolute_path)
                    logger.error("---------")
                    logger.error(e.absolute_schema_path)
                    exit(-1)
