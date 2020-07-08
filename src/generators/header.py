import os
from datetime import datetime
from string import Template
from common import get_ordered_species_dict
from common import ordered_taxon_species_map_from_data_dictionary


class HeaderTemplate(Template):

    delimiter = '%'


def create_header(file_type, database_version, data_format, config_info='', readme='', stringency_filter='', taxon_ids=''):
    """

    :param file_type:
    :param database_version:
    :param data_format:
    :param config_info: if not provide will pull species from YAML data dictionary file
    :param readme:
    :param stringency_filter:
    :param taxon_ids:
    :param file_type:
    :return:
    """

    ordered_taxon_species_map = get_ordered_species_dict(config_info, taxon_ids) if config_info != '' \
        else ordered_taxon_species_map_from_data_dictionary(taxon_ids)

    if stringency_filter != '':
        stringency_filter = '\n# Orthology Filter: ' + stringency_filter

    gen_time = datetime.utcnow().strftime("%Y-%m-%d %H:%M")
    metadata = {'filetype': file_type,
                'databaseVersion': database_version,
                'genTime': gen_time,
                'stringencyFilter': stringency_filter,
                'dataFormat': data_format,
                'readme': readme}

    if data_format == 'json':
        species = []
        for key, value in ordered_taxon_species_map.items():
            species.append({'taxonId': key,
                            'speciesName': value})
        metadata['species'] = species

        return metadata
    elif data_format == 'tsv':
        metadata['taxonIds'] = ', '.join(ordered_taxon_species_map.keys())
        metadata['species'] = ', '.join(ordered_taxon_species_map.values())

        my_path = os.path.abspath(os.path.dirname(__file__))
        file_header_template = HeaderTemplate(open(my_path + '/header_template.txt').read())

        return file_header_template.substitute(metadata)
    else:
        raise ValueError("Wrong data_format: " + "' - must be set to 'json' or 'tsv'")
