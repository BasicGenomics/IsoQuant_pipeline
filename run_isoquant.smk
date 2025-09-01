configfile: "config/config.yaml"

if 'resource_dir' not in config:
    config['resource_dir'] = '/usr/local/app/resources/'
else:
    if config['resource_dir'] == '':
        config['resource_dir'] = '/usr/local/app/resources/'


REFFILE = "{resource_dir}genome_references/{species}/reference.fa".format(resource_dir = config["resource_dir"], species = config["species"])
GTFFILE = "{resource_dir}genome_references/{species}/geneannotations".format(resource_dir=config["resource_dir"], species = config["species"])

rule all:
    input: "results/isoquant_output/mudata/{project_name}_counts.h5mu".format(project_name=config["project_name"])

rule modify_bamfile:
    input: "results/{name}.stitched.molecules.sorted.bam".format(project_name=config["project_name"])
    output: temp("results/isoquant_output/{project_name}.isoquant_modified.bam".format(project_name=config["project_name"]))
    shell: "python scripts/modify_file_for_IsoQuant.py --input {input} --output {output}"

rule index_bamfile:
    input: "results/isoquant_output/{project_name}.isoquant_modified.bam".format(project_name=config["project_name"])
    output:  temp("results/isoquant_output/{project_name}.isoquant_modified.bam.bai".format(project_name=config["project_name"]))
    shell: "samtools index {input}"

if config['bulk']:
    rule run_IsoQuant:
        input: bam =  "results/isoquant_output/{project_name}.isoquant_modified.bam".format(project_name = config["project_name"]), bai =  "results/isoquant_output/{project_name}.isoquant_modified.bam.bai".format(project_name = config["project_name"])
        output: counts = "results/isoquant_output/{project_name}/{project_name}.transcript_model_grouped_counts_linear.tsv".format(project_name = config["project_name"]), models = "results/isoquant_output/{project_name}/{project_name}.transcript_models.gtf".format(project_name = config["project_name"])
        params: fasta = REFFILE , gtf = "{}.gff3".format(GTFFILE)
        shell: "python IsoQuant/isoquant.py --reference {params.fasta} --genedb {params.gtf} --bam {input.bam} --data_type pacbio_ccs -p {config[project_name]} -o results/isoquant_output/ --read_group tag:SM --complete_genedb --count_exons"
else:
    rule run_IsoQuant:
        input: bam =  "results/isoquant_output/{project_name}.isoquant_modified.bam".format(project_name = config["project_name"]), bai =  "results/isoquant_output/{project_name}.isoquant_modified.bam.bai".format(project_name = config["project_name"])
        output: counts = "results/isoquant_output/{project_name}/{project_name}.transcript_model_grouped_counts_linear.tsv".format(project_name = config["project_name"]), models = "results/isoquant_output/{project_name}/{project_name}.transcript_models.gtf".format(project_name = config["project_name"])
        params: fasta = REFFILE , gtf = "{}.gff3".format(GTFFILE)
        shell: "python IsoQuant/isoquant.py --reference {params.fasta} --genedb {params.gtf} --bam {input.bam} --data_type pacbio_ccs -p {config[project_name]} -o results/isoquant_output/ --read_group tag:BC --complete_genedb --count_exons"

rule make_MuData:
    input: counts = "results/isoquant_output/{project_name}/{project_name}.transcript_model_grouped_counts_linear.tsv".format(project_name = config["project_name"]), models = "results/isoquant_output/{project_name}/{project_name}.transcript_models.gtf".format(project_name = config["project_name"])
    output: "results/isoquant_output/mudata/{project_name}_counts.h5mu".format(project_name=config["project_name"])
    shell: "python scripts/make_mudata.py --counts {input.counts} --models {input.models} --output {output}"
