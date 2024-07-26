import os
import importlib.resources
from stemcnv_check import STEM_CNV_CHECK

rule filter_snp_vcf:
  input: os.path.join(DATAPATH, "{sample_id}", "{sample_id}.unprocessed.vcf")
  output: pipe(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.processed-SNP-data.{filter}-filter.vcf"))
  threads: get_tool_resource('filter_snp_vcf', 'threads')
  resources:
    runtime=get_tool_resource('filter_snp_vcf', 'runtime'),
    mem_mb=get_tool_resource('filter_snp_vcf', 'memory'),
    partition=get_tool_resource('filter_snp_vcf', 'partition')
  log:
    err=os.path.join(LOGPATH, "filter_snp_vcf", "{sample_id}", "{filter}.error.log"),
    #out=os.path.join(LOGPATH, "filter_snp_vcf", "{sample_id}", "{filter}.out.log")
  conda:
    importlib.resources.files(STEM_CNV_CHECK).joinpath("envs","python-vcf.yaml")
  script:
    '../scripts/filter_snp_vcf.py'
    

rule annotate_snp_vcf:
    input:
      vcf = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.processed-SNP-data.{filter}-filter.vcf"),
      genomefasta =  config['static_data']['genome_fasta_file']
    output: os.path.join(DATAPATH, "{sample_id}", "{sample_id}.annotated-SNP-data.{filter}-filter.vcf.gz")
    threads: get_tool_resource('annotate_snp_vcf', 'threads')
    resources:
        runtime=get_tool_resource('annotate_snp_vcf', 'runtime'),
        mem_mb=get_tool_resource('annotate_snp_vcf', 'memory'),
        partition=get_tool_resource('annotate_snp_vcf', 'partition')
    log:
        err=os.path.join(LOGPATH, "annotate_snp_vcf", "{sample_id}", "{filter}.error.log"),
        #out=os.path.join(LOGPATH, "annotate_snp_vcf", "{sample_id}", "{filter}.out.log")
    params:
        genomeversion = 'GRCh38' if config['genome_build'] == 'hg38' else 'GRCh37',
        vep_cache_path = config['static_data']['VEP_cache_path'],
    conda:
      importlib.resources.files(STEM_CNV_CHECK).joinpath("envs", "vep-annotation.yaml")
    shell:
        'vep --verbose ' 
        '--fasta {input.genomefasta} ' # only needed for HGVS
        '--input_file {input.vcf} '
        '--output_file {output} '
        '--compress_output bgzip '
        '--vcf '
        '--force_overwrite '
        '--assembly {params.genomeversion} '

        #'--cache '
        #'--cache_version 112 '
        #'--dir_cache {params.vep_cache_path} '
        #> possibly useful to figure out protein coverage
        '--total_length '
        #'--numbers' # would add EXON,INTRON fields
        # Gene & protein annotation
        '--gencode_basic ' #> limit to gencode transcripts
        '--symbol ' #> add gene symbol
        '--terms SO ' # how to write/format/annotate the consequence
        '--hgvs ' #> add HGVS nomenclature (protrein changes)
        # maybe: '--no_esacpe' applies to HGVS
        '--pick ' #> pick the most severe consequence (& gene) per variant
        
        #> can check for existing annotations from ClinVar, COSMIC etc
        '--gene_phenotype '
        # will query databases for existing annotation at the same position
        '--check_existing --no_check_alleles '
        '--af ' # af only maybe useful?
        ##> probably overkill (enhancers/promotors etc 
        #'--regulatory'
        # Select content of CSQ, reduced need to drop stuff after
        # The first 3 are only for tetsing / confirmation
        # Probably remove: Feature,Feature_type
        # Removed: Amino_acids,Codons
        # How-to-type?: cDNA_position,CDS_position,Protein_position
        '--fields "Uploaded_variation,Location,Allele,Gene,SYMBOL,STRAND,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Extra,HGVSc,HGVSp,GENE_PHENO,Existing_variation,CLIN_SIG,SOMATIC,PHENO,AF" '
        # # Extra scores for proten changes (from SNPs)
        # #> probably not needed for array, unless some probes are on 'unknown' positions
        # '--polyphen b'
        # '--sift b'
        # # plugins: CADD, revel, AlphaMisssense
        ' 2> {log.err}'