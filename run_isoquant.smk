configfile: "config/config.yaml"


REFFILE = "{resource_dir}genome_references/{species}/current/reference.fa".format(resource_dir = config["resource_dir"], species = config["species"])
GTFFILE = "{resource_dir}genome_references/{species}/current/geneannotations".format(resource_dir=config["resource_dir"], species=config["species"])
rule all:
    input: "results/isoquant_output/mudata.done"

rule run_IsoQuant:
    input: "results/{project_name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.molecules.sorted.bam".format(project_name=config["project_name"])
    output: "results/isoquant_output/{project_name}.transcript_model_grouped_counts.tsv".format(project_name = config["project_name"]), "results/isoquant_output/{project_name}.transcript_models.gtf".format(project_name = config["project_name"])
    conda: "requirements.yaml"
    params: fasta = REFFILE , gtf = "{}.isoquant.gtf".format(GTFFILE)
    shell: "python {workflow.basedir}/IsoQuant/isoquant.py --reference {params.fasta} --genedb {params.gtf} --bam {input} --data_type pacbio_ccs -p {config[project_name]} -o results/isoquant_output/ --read_group tag:SM --complete_genedb --count_exons"

rule make_MuData:
    input: counts = "results/isoquant_output/{project_name}/{project_name}.transcript_model_grouped_counts.tsv".format(project_name = config["project_name"]), models = "results/isoquant_output/{project_name}/{project_name}.transcript_models.gtf".format(project_name = config["project_name"])
    output: touch("results/isoquant_output/mudata.done")
    conda: "requirements.yaml"
    shell: "python {workflow.basedir}/make_mudata.py --counts {input.counts} --models {input.models} --output-dir results/isoquant_output/mudata/"
