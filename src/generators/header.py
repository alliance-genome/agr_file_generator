
import os
from datetime import datetime
from string import Template

class HeaderTemplate(Template):

    delimiter = '%'


def create_header(file_type, database_version, stringency_filter='', taxon_ids=''):
    """

    :param file_type:
    :param database_version:
    :param stringency_filter:
    :param taxon_ids:
    :return:
    """

    if stringency_filter != '':
        stringency_filter = '# Orthology Filter: ' + stringency_filter

    gen_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    to_change = {'filetype': file_type, 'taxon_ids': taxon_ids, 'database_version': database_version,
                 'gen_time': gen_time, 'stringency_filter': stringency_filter}

    my_path = os.path.abspath(os.path.dirname(__file__))
    file_header_template = HeaderTemplate(open(my_path + '/template.txt').read())

    file_header = file_header_template.substitute(to_change)

    return file_header