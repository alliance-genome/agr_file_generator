import logging
import os

import wget


logger = logging.getLogger(name=__name__)


def get_local_path(local_assembly_dir, assembly):
    return os.path.join(local_assembly_dir, assembly + '.fa')


def download(local_assembly_dir, s3_uri, assembly):
    url = s3_uri + assembly + '.fa'
    os.makedirs(local_assembly_dir, exist_ok=True)
    download_to = get_local_path(local_assembly_dir, assembly)
    logger.info('Downloading assembly sequence (FASTA) from %s to %s',
                url,
                download_to)
    wget.download(url, out=local_assembly_dir)


def ensure_downloaded(local_assembly_dir, assembly_to_s3_uri):
    for (ass_name, s3_uri) in assembly_to_s3_uri.items():
        local_path = get_local_path(local_assembly_dir, ass_name)
        if os.path.exists(local_path):
            logger.warning('Assembly Sequence already downloaded for %r, skipping',
                           ass_name)
        else:
            download(local_assembly_dir, s3_uri, ass_name)


def get_seq(fasta, chromosome, start, end):
    start -= 1
    if chromosome.startswith('chr'):
        chromosome_str = chromosome[3:]
    else:
        chromosome_str = chromosome
    return fasta[chromosome_str][start:end].seq
