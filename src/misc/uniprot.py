# Paulo Nuin Feb 2020

import sys
import json
import csv
import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG')

if __name__ == '__main__':

    logger.info('Reading ' + sys.argv[1])
    json_input = json.load(open(sys.argv[1]))
    logger.info('Reading completed')
    uniprot_output = open('uniprot.tsv', 'w')
    uniprot_output.write('UniProtID\tIdentifier\n')
    uniprot_output.close()

    # for i in json_input:
    #     print(i['GlobalCrossReferenceID'].find('UniProt'))

    logger.info('Printing output file uniprot.tsv')
    with open('uniprot.tsv', 'a') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        for item in json_input:
            if item['GlobalCrossReferenceID'].find('UniProt') >= 0:
                writer.writerow([item['GlobalCrossReferenceID'].replace('UniProtKB:', ''), item['GeneID']])
    logger.info('Output completed')

