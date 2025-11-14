configfile: "config/config.yaml"

if 'resource_dir' not in config:
    config['resource_dir'] = '/usr/local/app/resources/'
else:
    if config['resource_dir'] == '':
        config['resource_dir'] = '/usr/local/app/resources/'


REFFILE = "{resource_dir}genome_references/{reference}/reference.fa".format(resource_dir = config["resource_dir"], reference = config["reference"])
GTFFILE = "{resource_dir}genome_references/{reference}/geneannotations".format(resource_dir=config["resource_dir"], reference = config["reference"])

rule all:
    input: "results/isoquant_output/mudata/{name}_counts.h5mu".format(name=config["name"])

rule modify_bamfile:
    input: "results/{name}.stitched.molecules.sorted.bam".format(name=config["name"])
    output: temp("results/isoquant_output/{name}.isoquant_modified.bam".format(name=config["name"]))
    shell: "python scripts/modify_file_for_IsoQuant.py --input {input} --output {output}"

rule index_bamfile:
    input: "results/isoquant_output/{name}.isoquant_modified.bam".format(name=config["name"])
    output:  temp("results/isoquant_output/{name}.isoquant_modified.bam.bai".format(name=config["name"]))
    shell: "samtools index {input}"

rule run_IsoQuant:
    input: bam =  "results/isoquant_output/{name}.isoquant_modified.bam".format(name = config["name"]), bai =  "results/isoquant_output/{name}.isoquant_modified.bam.bai".format(name = config["name"])
    output: counts = "results/isoquant_output/{name}/{name}.transcript_model_grouped_counts_linear.tsv".format(name = config["name"]), models = "results/isoquant_output/{name}/{name}.transcript_models.gtf".format(name = config["name"])
    params: fasta = REFFILE , gtf = "{}.gff3".format(GTFFILE)
    shell: "python IsoQuant/isoquant.py --reference {params.fasta} --genedb {params.gtf} --bam {input.bam} --data_type pacbio_ccs -p {config[name]} -o results/isoquant_output/ --read_group tag:SM --complete_genedb --count_exons"

rule make_MuData:
    input: counts = "results/isoquant_output/{name}/{name}.transcript_model_grouped_counts_linear.tsv".format(name = config["name"]), models = "results/isoquant_output/{name}/{name}.transcript_models.gtf".format(name = config["name"])
    output: "results/isoquant_output/mudata/{name}_counts.h5mu".format(name=config["name"])
    shell: "python scripts/make_mudata.py --counts {input.counts} --models {input.models} --output {output}"
