import os
from datetime import datetime
from string import Template


class HeaderTemplate(Template):

    delimiter = '%'


def create_header(file_type, database_version, data_format='', readme='', stringency_filter='', ordered_taxon_species_map=''):
    """

    :param file_type:
    :param database_version:
    :param stringency_filter:
    :param taxon_ids:
    :return:
    """

    if stringency_filter != '':
        stringency_filter = '\n# Orthology Filter: ' + stringency_filter

    gen_time = datetime.utcnow().strftime("%Y-%m-%d %H:%M")
    to_change = {"filetype": file_type,
                 "taxon_ids": ", ".join(ordered_taxon_species_map.keys()),
                 "species":  ", ".join(ordered_taxon_species_map.values()),
                 'database_version': database_version,
                 'gen_time': gen_time,
                 'stringency_filter': stringency_filter,
                 'data_format': data_format,
                 'readme': readme}

    my_path = os.path.abspath(os.path.dirname(__file__))
    file_header_template = HeaderTemplate(open(my_path + '/header_template.txt').read())

    file_header = file_header_template.substitute(to_change)

    return file_header
