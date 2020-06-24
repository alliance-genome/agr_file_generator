
import logging
import os
import time
import click
import coloredlogs
from common import ContextInfo
from common import get_neo_uri
from data_source import DataSource
from generators import (daf_file_generator, db_summary_file_generator,
                        expression_file_generator,
                        gene_cross_reference_file_generator,
                        orthology_file_generator, vcf_file_generator,
                        uniprot_cross_reference_generator)


config_info = ContextInfo()
debug_level = logging.DEBUG if config_info.config["DEBUG"] else logging.INFO
neo_debug_level = logging.DEBUG if config_info.config["NEO_DEBUG"] else logging.INFO

if config_info.config["GENERATED_FILES_FOLDER"]:
    generated_files_folder = config_info.config["GENERATED_FILES_FOLDER"]
else:
    generated_files_folder = os.path.join("/tmp", "agr_generated_files")

coloredlogs.install(level=debug_level,
                    fmt='%(asctime)s %(levelname)s: %(name)s:%(lineno)d: %(message)s',
                    field_styles={'asctime': {'color': 'green'},
                                  'hostname': {'color': 'magenta'},
                                  'levelname': {'color': 'white', 'bold': True},
                                  'name': {'color': 'blue'},
                                  'programname': {'color': 'cyan'}})


logging.getLogger("urllib3").setLevel(debug_level)
logging.getLogger("neobolt").setLevel(neo_debug_level)
logger = logging.getLogger(__name__)

if config_info.config["DEBUG"]:
    logger.warning('DEBUG mode enabled!')

taxon_id_fms_subtype_map = {"NCBITaxon:10116": "RGD",
                            "NCBITaxon:9606": "HUMAN",
                            "NCBITaxon:7227": "FB",
                            "NCBITaxon:6239": "WB",
                            "NCBITaxon:7955": "ZFIN",
                            "NCBITaxon:10090": "MGI",
                            "NCBITaxon:559292": "SGD"}


@click.command()
@click.option('--vcf', is_flag=True, help='Generates VCF files')
@click.option('--orthology', is_flag=True, help='Generates orthology files')
@click.option('--disease', is_flag=True, help='Generates DAF files')
@click.option('--expression', is_flag=True, help='Generates expression files')
@click.option('--db-summary', is_flag=True, help='Generates summary of database contents')
@click.option('--gene-cross-reference', is_flag=True, help='Generates a file of cross references for gene objects')
@click.option('--all-filetypes', is_flag=True, help='Generates all filetypes')
@click.option('--tab', is_flag=True, help='Generates tab delimited files with VCF info columns contents')
@click.option('--uniprot', is_flag=True, help='Generates a file of gene and uniprot cross references')
@click.option('--upload', is_flag=True, help='Submits generated files to File Management System (FMS)')
def main(vcf,
         orthology,
         disease,
         expression,
         db_summary,
         gene_cross_reference,
         all_filetypes,
         upload,
         tab,
         uniprot,
         generated_files_folder=os.path.abspath(os.path.join(os.getcwd(), os.pardir)) + '/output',
         input_folder=os.path.abspath(os.path.join(os.getcwd(), os.pardir)) + '/input',
         fasta_sequences_folder='sequences',
         skip_chromosomes={'Unmapped_Scaffold_8_D1580_D1567'}):

    start_time = time.time()
    click.echo("Start Time: %s" % time.strftime("%H:%M:%S", time.gmtime(start_time)))

    if not os.path.exists(generated_files_folder):
        os.makedirs(generated_files_folder)

    click.echo('INFO:\tFiles output: ' + generated_files_folder)
    if vcf is True or all_filetypes is True:
        if not tab:
            click.echo('INFO:\tGenerating VCF files')
        else:
            click.echo('INFO:\tGenerating TAB delimited files')
        generate_vcf_files(generated_files_folder, fasta_sequences_folder, skip_chromosomes, config_info, upload, tab)
    if orthology is True or all_filetypes is True:
        click.echo('INFO:\tGenerating Orthology file')
        generate_orthology_file(generated_files_folder, config_info, upload)
    if disease is True or all_filetypes is True:
        click.echo('INFO:\tGenerating Disease files')
        generate_daf_file(generated_files_folder, config_info, taxon_id_fms_subtype_map, upload)
    if expression is True or all_filetypes is True:
        click.echo('INFO:\tGenerating Expression files')
        generate_expression_file(generated_files_folder, config_info, taxon_id_fms_subtype_map, upload)
    if db_summary is True or all_filetypes is True:
        click.echo('INFO:\tGenerating DB summary file')
        generate_db_summary_file(generated_files_folder, config_info, upload)
    if gene_cross_reference is True or all_filetypes is True:
        click.echo('INFO:\tGenerating Gene Cross Reference file')
        generate_gene_cross_reference_file(generated_files_folder, config_info, upload)
    if uniprot is True or all_filetypes is True:
        generate_uniprot_cross_reference(generated_files_folder, input_folder, config_info, upload)

    end_time = time.time()
    elapsed_time = end_time - start_time
    click.echo('File Generator finished. Elapsed time: %s' % time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))


def generate_vcf_file(assembly, generated_files_folder, fasta_sequence_folder, skip_chromosomes, config_info, upload_flag, tab_flag):
    logger.info("Querying Assembly: " + assembly)

    variants_query = '''MATCH (s:Species)-[:FROM_SPECIES]-(a:Allele)-[:VARIATION]-(v:Variant)-[l:LOCATED_ON]->(c:Chromosome),
                              (v:Variant)-[:VARIATION_TYPE]->(st:SOTerm),
                              (v:Variant)-[:ASSOCIATION]->(p:GenomicLocation)-[:ASSOCIATION]->(assembly:Assembly {primaryKey: "''' + assembly + '''"})
                     WHERE NOT v.genomicReferenceSequence = v.genomicVariantSequence
                           OR v.genomicVariantSequence = ""
                     OPTIONAL MATCH (a:Allele)-[:IS_ALLELE_OF]-(g:Gene)
                     WITH COLLECT(DISTINCT {symbol: a.symbol,
                                            symbolText: a.symbolText,
                                            id: a.primaryKey}) AS alleles,
                          s, v, l, c, st, p, assembly
                     OPTIONAL MATCH (v:Variant)-[:ASSOCIATION]-(glc:GeneLevelConsequence)-[:ASSOCIATION]-(g:Gene)
                     OPTIONAL MATCH (v:Variant)-[:ASSOCIATION]-(tlc:TranscriptLevelConsequence)-[:ASSOCIATION]-(t:Transcript)
                     RETURN c.primaryKey AS chromosome,
                            v.globalId AS globalId,
                            right(v.paddingLeft,1) AS paddingLeft,
                            v.genomicReferenceSequence AS genomicReferenceSequence,
                            v.genomicVariantSequence AS genomicVariantSequence,
                            v.hgvsNomenclature AS hgvsNomenclature,
                            v.dataProvider AS dataProvider,
                            assembly.primaryKey AS assembly,
                            alleles,
                            COLLECT(DISTINCT {gene: g.primaryKey,
                                              geneSymbol: g.symbol,
                                              consequence: glc.geneLevelConsequence,
                                              impact: glc.impact}) AS geneConsequences,
                            collect(DISTINCT {transcript: t.primaryKey,
                                              transcriptGFF3ID: t.gff3ID,
                                              transcriptGFF3Name: t.name,
                                              consequence: tlc.transcriptLevelConsequence,
                                              impact: tlc.impact}) AS transcriptConsequences,
                            p.start AS start,
                            p.end AS end,
                            s.name AS species,
                            st.nameKey AS soTerm
                     '''

    if config_info.config["DEBUG"]:
        logger.info(variants_query)
        start_time = time.time()
        logger.info("Start time: %s", time.strftime("%H:%M:%S", time.gmtime(start_time)))

    data_source = DataSource(get_neo_uri(config_info), variants_query)
    gvf = vcf_file_generator.VcfFileGenerator(data_source,
                                              generated_files_folder,
                                              config_info)
    gvf.generate_files(skip_chromosomes=skip_chromosomes, upload_flag=upload_flag, tab_flag=tab_flag)

    if config_info.config["DEBUG"]:
        end_time = time.time()
        logger.info("Created VCF file - End time: %s", time.strftime("%H:%M:%S", time.gmtime(end_time)))
        logger.info("Time Elapsed: %s", time.strftime("%H:%M:%S", time.gmtime(end_time - start_time)))


def generate_vcf_files(generated_files_folder, fasta_sequences_folder, skip_chromosomes, config_info, upload_flag, tab_flag):
    os.makedirs(generated_files_folder, exist_ok=True)
    os.makedirs(fasta_sequences_folder, exist_ok=True)

    assembly_query = """MATCH (a:Assembly)
                        RETURN a.primaryKey as assemblyID"""
    assembly_data_source = DataSource(get_neo_uri(config_info), assembly_query)

    if config_info.config["DEBUG"]:
        start_time = time.time()
        logger.info("Start time for generating VCF files: %s", time.strftime("%H:%M:%S", time.gmtime(start_time)))

    for assembly_result in assembly_data_source:
        assembly = assembly_result["assemblyID"]
        if assembly not in ["", "GRCh38", "R64-2-1"]:
            generate_vcf_file(assembly,
                              generated_files_folder,
                              fasta_sequences_folder,
                              skip_chromosomes,
                              config_info,
                              upload_flag,
                              tab_flag)

    if config_info.config["DEBUG"]:
        end_time = time.time()
        logger.info("Created VCF files - End time: %s", time.strftime("%H:%M:%S", time.gmtime(end_time)))
        logger.info("Time Elapsed: %s", time.strftime("%H:%M:%S", time.gmtime(end_time - start_time)))


def generate_orthology_file(generated_files_folder, config_info, upload_flag):
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

    if config_info.config["DEBUG"]:
        logger.info("Orthology query")
        logger.info(orthology_query)
        start_time = time.time()
        logger.info("Start time: %s", time.strftime("%H:%M:%S", time.gmtime(start_time)))

    data_source = DataSource(get_neo_uri(config_info), orthology_query)
    of = orthology_file_generator.OrthologyFileGenerator(data_source,
                                                         generated_files_folder,
                                                         config_info)
    of.generate_file(upload_flag=upload_flag)

    if config_info.config["DEBUG"]:
        end_time = time.time()
        logger.info("Created VCF file - End time: %s", time.strftime("%H:%M:%S", time.gmtime(end_time)))
        logger.info("Time Elapsed: %s", time.strftime("%H:%M:%S", time.gmtime(end_time - start_time)))


def generate_daf_file(generated_files_folder, config_info, taxon_id_fms_subtype_map, upload_flag):
    daf_query = '''MATCH (disease:DOTerm)-[:ASSOCIATION]-(dej:Association:DiseaseEntityJoin)-[:ASSOCIATION]-(object)-[:FROM_SPECIES]-(species:Species)
                   WHERE (object:Gene OR object:Allele OR object:AffectedGenomicModel)
                         AND dej.joinType IN ["IS_MARKER_FOR", // need to remove when removed from database
                                              "IS_IMPLICATED_IN", // need to remove when removed from database
                                              "IS_MODEL_OF",
                                              "is_model_of",
                                              "is_implicated_in",
                                              "is_biomarker_for",
                                              "implicated_via_orthology",
                                              "biomarker_via_orthology"]
                   MATCH (dej:Association:DiseaseEntityJoin)-[:EVIDENCE]->(pj:PublicationJoin),
                         (p:Publication)-[:ASSOCIATION]->(pj:PublicationJoin)-[:ASSOCIATION]->(ec:Ontology:ECOTerm)
                   OPTIONAL MATCH (object:Gene)-[:ASSOCIATION]->(dej:Association:DiseaseEntityJoin)<-[:ASSOCIATION]-(otherAssociatedEntity)
                   OPTIONAL MATCH (pj:PublicationJoin)-[:MODEL_COMPONENT|PRIMARY_GENETIC_ENTITY]-(inferredFromEntity)
                   OPTIONAL MATCH (dej:Association:DiseaseEntityJoin)-[:FROM_ORTHOLOGOUS_GENE]->(oGene:Gene),
                                  (gene:Gene)-[o:ORTHOLOGOUS]->(oGene:Gene)
                   WHERE o.strictFilter AND ec.primaryKey IN ["ECO:0000250", "ECO:0000266", "ECO:0000501"] // ISS, ISO, and IEA respectively
                   //OPTIONAL MATCH (object)-[IS_ALLELE_OF]->(gene:Gene)
                   RETURN DISTINCT
                          dej.primaryKey as dejID,
                          species.primaryKey AS taxonId,
                          species.name AS speciesName,
                          collect(DISTINCT oGene.primaryKey) AS withOrthologs,
                          labels(object) AS objectType,
                          object.primaryKey AS dbObjectID,
                          object.symbol AS dbObjectSymbol,
                          object.name AS dbObjectName,
                          toLower(dej.joinType) AS associationType,
                          //collect(DISTINCT gene.primaryKey) AS inferredGeneAssociation,
                          disease.doId AS DOID,
                          disease.name as DOtermName,
                          collect(DISTINCT {pubModID: p.pubModId,
                                            pubMedID: p.pubMedId,
                                            evidenceCode:ec.primaryKey,
                                            evidenceCodeName: ec.name,
                                            inferredFromEntity: inferredFromEntity,
                                            otherAssociatedEntityID: otherAssociatedEntity.primaryKey}) as evidence,
                          REDUCE(t = "1900-01-01", c IN collect(left(pj.dateAssigned, 10)) | CASE WHEN c > t THEN c ELSE t END) AS dateAssigned,
                          ///takes most recent date
                          dej.dataProvider AS dataProvider'''

    if config_info.config["DEBUG"]:
        logger.info("Disease Association Query: ")
        logger.info(daf_query)
        start_time = time.time()
        logger.info("Start time: %s", time.strftime("%H:%M:%S", time.gmtime(start_time)))

    data_source = DataSource(get_neo_uri(config_info), daf_query)
    daf = daf_file_generator.DafFileGenerator(data_source,
                                              generated_files_folder,
                                              config_info,
                                              taxon_id_fms_subtype_map)
    daf.generate_file(upload_flag=upload_flag)

    if config_info.config["DEBUG"]:
        end_time = time.time()
        logger.info("Created Disease Association file - End time: %s", time.strftime("%H:%M:%S", time.gmtime(end_time)))
        logger.info("Time Elapsed: %s", time.strftime("%H:%M:%S", time.gmtime(end_time - start_time)))


def generate_expression_file(generated_files_folder, config_info, taxon_id_fms_subtype_map, upload_flag):
    expression_query = '''MATCH (speciesObj:Species)<-[:FROM_SPECIES]-(geneObj:Gene)-[:ASSOCIATION]->(begej:BioEntityGeneExpressionJoin)--(term)
                          WITH {primaryKey: speciesObj.primaryKey, name: speciesObj.name} AS species,
                               {primaryKey: geneObj.primaryKey, symbol: geneObj.symbol, dataProvider: geneObj.dataProvider} AS gene,
                               begej,
                               COLLECT(term) AS terms
                          MATCH (begej:BioEntityGeneExpressionJoin)<-[:ASSOCIATION]-(exp:ExpressionBioEntity)-[a:ANATOMICAL_STRUCTURE|CELLULAR_COMPONENT|ANATOMICAL_SUB_SUBSTRUCTURE|CELLULAR_COMPONENT_QUALIFIER|ANATOMICAL_SUB_STRUCTURE_QUALIFIER|ANATOMICAL_STRUCTURE_QUALIFIER]->(ontology:Ontology)
                          RETURN species,
                                 gene,
                                 terms,
                                 begej.primaryKey as begejId,
                                 exp.whereExpressedStatement AS location,
                                 COLLECT({edge: type(a),
                                          primaryKey: ontology.primaryKey,
                                          name: ontology.name}) AS ontologyPaths'''

    if config_info.config["DEBUG"]:
        logger.info("Expression query")
        logger.info(expression_query)
        start_time = time.time()
        logger.info("Start time: %s", time.strftime("%H:%M:%S", time.gmtime(start_time)))

    data_source = DataSource(get_neo_uri(config_info), expression_query)
    expression = expression_file_generator.ExpressionFileGenerator(data_source,
                                                                   generated_files_folder,
                                                                   config_info,
                                                                   taxon_id_fms_subtype_map)
    expression.generate_file(upload_flag=upload_flag)

    if config_info.config["DEBUG"]:
        end_time = time.time()
        logger.info("Created Expression file - End time: %s", time.strftime("%H:%M:%S", time.gmtime(end_time)))
        logger.info("Time Elapsed: %s", time.strftime("%H:%M:%S", time.gmtime(end_time - start_time)))


def generate_db_summary_file(generated_files_folder, config_info, upload_flag):
    db_summary_query = '''MATCH (entity)
                          WITH labels(entity) AS entityTypes
                          RETURN count(entityTypes) AS frequency,
                          entityTypes'''

    if config_info.config["DEBUG"]:
        logger.info("DB Summary Query")
        logger.info(db_summary_query)
        start_time = time.time()
        logger.info("Start time: %s", time.strftime("%H:%M:%S", time.gmtime(start_time)))

    data_source = DataSource(get_neo_uri(config_info), db_summary_query)
    db_summary = db_summary_file_generator.DbSummaryFileGenerator(data_source,
                                                                  generated_files_folder,
                                                                  config_info)
    db_summary.generate_file(upload_flag=upload_flag)

    if config_info.config["DEBUG"]:
        end_time = time.time()
        logger.info("Created DB Summary file - End time: %s", time.strftime("%H:%M:%S", time.gmtime(end_time)))
        logger.info("Time Elapsed: %s", time.strftime("%H:%M:%S", time.gmtime(end_time - start_time)))


def generate_gene_cross_reference_file(generated_files_folder, config_info, upload_flag):
    gene_cross_reference_query = '''MATCH (g:Gene)--(cr:CrossReference)
                          RETURN g.primaryKey as GeneID,
                                 cr.globalCrossRefId as GlobalCrossReferenceID,
                                 cr.crossRefCompleteUrl as CrossReferenceCompleteURL,
                                 cr.page as ResourceDescriptorPage,
                                 g.taxonId as TaxonID'''

    if config_info.config["DEBUG"]:
        logger.info("Gene Cross Reference query")
        logger.info(gene_cross_reference_query)
        start_time = time.time()
        logger.info("Start time: %s", time.strftime("%H:%M:%S", time.gmtime(start_time)))

    data_source = DataSource(get_neo_uri(config_info), gene_cross_reference_query)
    gene_cross_reference = gene_cross_reference_file_generator.GeneCrossReferenceFileGenerator(data_source,
                                                                                               generated_files_folder,
                                                                                               config_info)
    gene_cross_reference.generate_file(upload_flag=upload_flag)

    if config_info.config["DEBUG"]:
        end_time = time.time()
        logger.info("Gene Cross Reference file - End time: %s", time.strftime("%H:%M:%S", time.gmtime(end_time)))
        logger.info("Time Elapsed: %s", time.strftime("%H:%M:%S", time.gmtime(end_time - start_time)))


def generate_uniprot_cross_reference(generated_files_folder, input_folder, config_info, upload_flag):
    uniprot_cross_reference_query = '''MATCH (g:Gene)--(cr:CrossReference)
                                WHERE cr.prefix = "UniProtKB"
                                RETURN g.primaryKey as GeneID,
                                    cr.globalCrossRefId as GlobalCrossReferenceID'''

    if config_info.config["DEBUG"]:
        logger.info("UniProt Cross Reference query")
        logger.info(uniprot_cross_reference_query)
        start_time = time.time()
        logger.info("Start time: %s", time.strftime("%H:%M:%S", time.gmtime(start_time)))

    data_source = DataSource(get_neo_uri(config_info), uniprot_cross_reference_query)
    ucf = uniprot_cross_reference_generator.UniProtGenerator(data_source, config_info, generated_files_folder)
    ucf.generate_file(upload_flag=upload_flag)

    if config_info.config["DEBUG"]:
        end_time = time.time()
        logger.info("Created UniProt Cross Reference file - End time: %s", time.strftime("%H:%M:%S", time.gmtime(end_time)))
        logger.info("Time Elapsed: %s", time.strftime("%H:%M:%S", time.gmtime(end_time - start_time)))


if __name__ == '__main__':
    main()
