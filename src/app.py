import logging
import os

import click

import VcfFileGenerator
import OrthologyFileGenerator
import DafFileGenerator
import ExpressionFileGenerator
from data_source import DataSource


host = os.environ.get('NEO4J_HOST', 'build.alliancegenome.org')
port = int(os.environ.get('NEO4J_PORT', 7687))
alliance_db_version = os.environ.get('ALLIANCE_RELEASE', '2.3.0')

uri = "bolt://" + host + ":" + str(port)


def setup_logging(logger_name):
    logging.basicConfig(level=logging.DEBUG)

# Common configuration variables used throughout the script.
class ContextInfo(object):

    def __init__(self):
        config_file = open('src/config.yaml', 'r')
        self.config = yaml.load(config_file, Loader=yaml.FullLoader)

        # Look for ENV variables to replace default variables from config file.
        for key in self.config.keys():
            try: 
                self.config[key] = os.environ[key]
            except KeyError:
                logger.info('Environmental variable not found for \'{}\'. Using config.yaml value.'.format(key))
                pass  # If we don't find an ENV variable, keep the value from the config file.
        
        logger.debug('Initialized with config values: {}'.format(self.config))
        logger.info('Retrieval errors will be emailed to: {}'.format(self.config['notification_emails'])) 

@click.command()
@click.option('--vcf', is_flag=True, help='Generates VCF files')
@click.option('--ortho', is_flag=True, help='Generates orthology files')
@click.option('--daf', is_flag=True, help='Generates DAF files')
@click.option('--expr', is_flag=True, help='Generates expression files')
@click.option('--upload', is_flag=True, help='Submits generated files to File Management System (FMS)')
def main(vcf, ortho, daf, expr, upload,
         generated_files_folder=os.path.abspath(os.path.join(os.getcwd(), os.pardir)) + '/output',
         fasta_sequences_folder='sequences',
         skip_chromosomes={'Unmapped_Scaffold_8_D1580_D1567'}):
    
    click.echo('INFO:\tFiles output: ' + generated_files_folder)
    context_info = ContextInfo()
    if vcf is True:
        click.echo('INFO:\tGenerating VCF files')
        generate_vcf_files(generated_files_folder, fasta_sequences_folder, skip_chromosomes, context_info, upload)
    elif ortho is True:
        click.echo('INFO:\tGenerating Orthology file')
        generate_orthology_file(generated_files_folder, alliance_db_version)
    elif daf is True:
        click.echo('INFO:\tGenerating DAF file')
        generate_daf_file(generated_files_folder, alliance_db_version)
    elif expr is True:
        click.echo('INFO:\tGenerating Expression file')
        generate_expression_file(generated_files_folder, alliance_db_version)

def generate_vcf_files(generated_files_folder, fasta_sequences_folder, skip_chromosomes, context_info, upload_flag):
    os.makedirs(generated_files_folder, exist_ok=True)
    os.makedirs(fasta_sequences_folder, exist_ok=True)
    variants_query = """MATCH (s:Species)-[:FROM_SPECIES]-(a:Allele)-[:VARIATION]-(v:Variant)-[l:LOCATED_ON]-(c:Chromosome),
                              (v:Variant)-[:VARIATION_TYPE]-(st:SOTerm),
                              (v:Variant)-[:ASSOCIATION]-(p:GenomicLocation)
                     OPTIONAL MATCH (a:Allele)-[:IS_ALLELE_OF]-(g:Gene)
                     RETURN c.primaryKey AS chromosome,
                            v.globalId AS globalId,
                            right(v.paddingLeft,1) AS paddingLeft,
                            v.genomicReferenceSequence AS genomicReferenceSequence,
                            v.genomicVariantSequence AS genomicVariantSequence,
                            v.geneLevelConsequence AS geneLevelConsequence,
                            v.hgvsNomenclature AS hgvsNomenclature,
                            v.dataProvider AS dataProvider,
                            a.symbol AS symbol,
                            a.symbolText as symbolText,
                            p.assembly AS assembly,
                            collect(a.primaryKey) AS alleles,
                            collect(g.primaryKey) AS geneSymbol,
                            CASE WHEN g IS NOT NULL THEN collect(g.primaryKey) ELSE [] END AS alleleOfGenes,
                            l.start AS start,
                            l.end AS end,
                            s.name AS species,
                            st.nameKey AS soTerm
                     """
    data_source = DataSource(uri, variants_query)
    gvf = VcfFileGenerator(data_source,
                           generated_files_folder,
                           context_info)
    gvf.generate_files(skip_chromosomes=skip_chromosomes, upload_flag=upload_flag)


def generate_orthology_file(generated_files_folder, alliance_db_version):
    orthology_query = '''MATCH (species1)<-[sa:FROM_SPECIES]-(gene1:Gene)-[o:ORTHOLOGOUS]->(gene2:Gene)-[sa2:FROM_SPECIES]->(species2:Species)
                       WHERE o.strictFilter
                       OPTIONAL MATCH (algorithm:OrthoAlgorithm)-[m:MATCHED]-(ogj:OrthologyGeneJoin)-[association:ASSOCIATION]-(gene1)
                       WHERE ogj.primaryKey = o.primaryKey
                       OPTIONAL MATCH (algorithm2:OrthoAlgorithm)-[m2:NOT_MATCHED]-(ogj2:OrthologyGeneJoin)-[ASSOCIATION]-(gene1)
                       WHERE ogj2.primaryKey = o.primaryKey
                       RETURN gene1.primaryKey AS gene1ID,
                              gene1.symbol AS gene1Symbol,
                              gene2.primaryKey AS gene2ID,
                              gene2.symbol AS gene2Symbol,
                              collect(DISTINCT algorithm.name) as Algorithms,
                              count(DISTINCT algorithm.name) AS numAlgorithmMatch,
                              count(DISTINCT algorithm2.name) AS numAlgorithmNotMatched,
                              toString(o.isBestScore) AS best,
                              toString(o.isBestRevScore) AS bestRev,
                              species1.primaryKey AS species1TaxonID,
                              species1.name AS species1Name,
                              species2.primaryKey AS species2TaxonID,
                              species2.name AS species2Name'''
    data_source = DataSource(uri, orthology_query)
    of = OrthologyFileGenerator(data_source,
                                generated_files_folder,
                                alliance_db_version)
    of.generate_file()


def generate_daf_file(generated_files_folder, alliance_db_version):
    daf_query = '''MATCH (dej:Association:DiseaseEntityJoin)-[:ASSOCIATION]-(object)-[da:IS_MARKER_FOR|:IS_IMPLICATED_IN|:IMPLICATED_VIA_ORTHOLOGY|:BIOMARKER_VIA_ORTHOLOGY]->(disease:DOTerm)
                   WHERE (object:Gene OR object:Allele)
                   AND da.uuid = dej.primaryKey
                   MATCH (object)-[FROM_SPECIES]->(species:Species)
                   OPTIONAL MATCH (ec:Ontology:ECOTerm)-[:ASSOCIATION]-(:PublicationEvidenceCodeJoin)-[:EVIDENCE]-(dej:Association:DiseaseEntityJoin)
                   OPTIONAL MATCH (p:Publication)-[:ASSOCIATION]-(:PublicationEvidenceCodeJoin)-[:EVIDENCE]-(dej:Association:DiseaseEntityJoin)
                   OPTIONAL MATCH (object)-[o:ORTHOLOGOUS]-(oGene:Gene)
                   WHERE o.strictFilter AND (ec.primaryKey = "ECO:0000250" OR ec.primaryKey = "ECO:0000266") // ISS and ISO respectively
                   OPTIONAL MATCH (object)-[IS_ALLELE_OF]->(gene:Gene)
                   RETURN  object.taxonId AS taxonId,
                           species.name AS speciesName,
                           collect(DISTINCT oGene.primaryKey) AS withOrthologs,
                           labels(object) AS objectType,
                           object.primaryKey AS dbObjectID,
                           object.symbol AS dbObjectSymbol,
                           p.pubMedId AS pubMedID,
                           p.pubModId As pubModID,
                           type(da) AS associationType,
                           collect(DISTINCT gene.primaryKey) AS inferredGeneAssociation,
                           disease.doId AS DOID,
                           disease.name as DOname,
                           ec.primaryKey AS evidenceCode,
                           dej.dateAssigned AS dateAssigned,
                           da.dataProvider AS dataProvider'''
    data_source = DataSource(uri, daf_query)
    daf = DafFileGenerator(data_source,
                           generated_files_folder,
                           alliance_db_version)
    daf.generate_file() 


def generate_expression_file(generated_files_folder, alliance_db_version):
    expression_query = '''MATCH (speciesObj:Species)<-[:FROM_SPECIES]-(geneObj:Gene)-[:ASSOCIATION]->(begej:BioEntityGeneExpressionJoin)--(term)
                          WITH {primaryKey: speciesObj.primaryKey, name: speciesObj.name} AS species,
                               {primaryKey: geneObj.primaryKey, symbol: geneObj.symbol} AS  gene,
                               begej,
                               COLLECT(term) AS terms
                          MATCH (begej:BioEntityGeneExpressionJoin)<-[:ASSOCIATION]-(exp:ExpressionBioEntity)-[a:ANATOMICAL_STRUCTURE|CELLULAR_COMPONENT|ANATOMICAL_SUB_SUBSTRUCTURE|CELLULAR_COMPONENT_QUALIFIER|ANATOMICAL_SUB_STRUCTURE_QUALIFIER|ANATOMICAL_STRUCTURE_QUALIFIER]->(ontology:Ontology)
                          //WHERE gene.primaryKey = 'ZFIN:ZDB-GENE-110411-206'
                          RETURN species, gene, terms, begej.primaryKey as begejId, exp.whereExpressedStatement as location,
                                                       COLLECT({edge: type(a),
                                                                primaryKey: ontology.primaryKey,
                                                                name: ontology.name}) as ontologyPaths'''
    data_source = DataSource(uri, expression_query)
    expression = ExpressionFileGenerator(data_source,
                                         generated_files_folder,
                                         alliance_db_version)
    expression.generate_file()


if __name__ == '__main__':
    main()
