from datetime import datetime
from string import Template

class HeaderTemplate(Template):
    """
    TBA
    """

    delimiter = '%'


def create_header(file_type, taxon_ids, database_version):

    gen_time = datetime.now()

    to_change = {'filetype': file_type, 'taxon_ids': taxon_ids, 'database_version': database_version, 'gen_time':gen_time}

    file_header_template = HeaderTemplate("""#########################################################################
#                                                                   
# %filetype                                                         
# Source: Alliance of Genome Resources (AGR)                        
# Orthology Filter: Stringent                                       
# TaxonIDs: %taxon_ids                                              
# Datebase Version: %database_version                               
# Date: %gen_time                                                   
#                                                                   
#########################################################################
""")

    file_header = file_header_template.substitute(to_change)

    print(file_header)

    return file_header