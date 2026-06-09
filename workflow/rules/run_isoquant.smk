rule modify_bamfile:
    input: "results/{name}.stitched.molecules.sorted.bam".format(name=config["name"])
    output: temp("results/isoquant/{name}.isoquant.bam".format(name=config["name"]))
    log: "results/isoquant/logs/{name}.modify_bamfile.log".format(name=config["name"])
    shell: "python workflow/scripts/modify_file_for_IsoQuant.py --input {input} --output {output} > {log} 2>&1"

rule index_bamfile:
    input: "results/isoquant/{name}.isoquant.bam".format(name=config["name"])
    output: temp("results/isoquant/{name}.isoquant.bam.bai".format(name=config["name"]))
    log: "results/isoquant/logs/{name}.index_bamfile.log".format(name=config["name"])
    shell: "samtools index {input} > {log} 2>&1"

rule run_isoquant:
    input: bam =  "results/isoquant/{name}.isoquant.bam".format(name = config["name"]), 
           bai =  "results/isoquant/{name}.isoquant.bam.bai".format(name = config["name"])
    output: counts = "results/isoquant/{name}/{name}.discovered_transcript_grouped_{token}_counts.linear.tsv".format(name = config["name"], token = GROUP_TOKEN),
            models = "results/isoquant/{name}/{name}.transcript_models.gtf".format(name = config["name"])
    params: ref = REF,
            gff3 = "{}.gff3".format(GFF),
    log: "results/isoquant/logs/{name}.run_isoquant.log".format(name=config["name"])
    shell: "python IsoQuant/isoquant.py --reference {params.ref} --genedb {params.gff3} --bam {input.bam} --data_type {config[data_type]} --basecode --basecode_max_gap {config[basecode_max_gap]} --matching_strategy {config[matching_strategy]} --stranded {config[stranded]} -p {config[name]} -o results/isoquant/ --read_group {config[read_group]} --complete_genedb --model_construction_strategy {config[model_construction_strategy]} --count_exons --transcript_quantification {config[transcript_quantification]} --gene_quantification {config[gene_quantification]} --polya_requirement {config[polya_requirement]} > {log} 2>&1"

rule make_mudata:
    input: counts = "results/isoquant/{name}/{name}.discovered_transcript_grouped_{token}_counts.linear.tsv".format(name = config["name"], token = GROUP_TOKEN),
           models = "results/isoquant/{name}/{name}.transcript_models.gtf".format(name = config["name"])
    output: "results/isoquant/mudata/{name}_counts.h5mu".format(name=config["name"])
    log: "results/isoquant/logs/{name}.make_mudata.log".format(name=config["name"])
    shell: "python workflow/scripts/make_mudata.py --counts {input.counts} --models {input.models} --output {output} > {log} 2>&1"