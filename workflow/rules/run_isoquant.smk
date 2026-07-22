rule modify_bam:
    input: "results/{name}.stitched.molecules.sorted.bam".format(name=config["name"])
    output: bam = temp("results/isoquant/{name}.isoquant.bam".format(name=config["name"])),
            tracking = "results/isoquant/{name}.adapted_molecules.tsv".format(name=config["name"])
    params: fl_flag = "--full-length-only" if config["full_length_only"] else ""
    log: "results/isoquant/logs/{name}.modify_bamfile.log".format(name=config["name"])
    shell: "python workflow/scripts/modify_bam.py --input {input} --output {output.bam} --tracking-file {output.tracking} {params.fl_flag} > {log} 2>&1"

rule index_bam:
    input: "results/isoquant/{name}.isoquant.bam".format(name=config["name"])
    output: temp("results/isoquant/{name}.isoquant.bam.bai".format(name=config["name"]))
    log: "results/isoquant/logs/{name}.index_bamfile.log".format(name=config["name"])
    shell: "samtools index {input} > {log} 2>&1"

rule collapse_annotation:
    input: "{}.gff3".format(GFF),
    output: COLLAPSED,
    log: "results/isoquant/logs/{name}.collapse_annotation.log".format(name=config["name"])
    shell: "python workflow/scripts/collapse_annotation.py --input {input} --output {output} --group-by {config[collapse_mode]} --representative {config[collapse_representative]} > {log} 2>&1"

def annotation_gff3(wildcards):
    return ANNOTATION

rule run_isoquant:
    input: bam =  "results/isoquant/{name}.isoquant.bam".format(name=config["name"]),
           bai =  "results/isoquant/{name}.isoquant.bam.bai".format(name=config["name"]),
           gff3 = annotation_gff3
    output: disc_transcript_counts = "results/isoquant/{name}/{name}.discovered_transcript_grouped_{token}_counts.linear.tsv".format(name=config["name"], token = GROUP_TOKEN),
            disc_transcript_tpm = "results/isoquant/{name}/{name}.discovered_transcript_grouped_{token}_tpm.tsv".format(name=config["name"], token = GROUP_TOKEN),
            disc_gene_counts = "results/isoquant/{name}/{name}.discovered_gene_grouped_{token}_counts.tsv".format(name=config["name"], token = GROUP_TOKEN),
            disc_gene_tpm = "results/isoquant/{name}/{name}.discovered_gene_grouped_{token}_tpm.tsv".format(name=config["name"], token = GROUP_TOKEN),
            ref_transcript_counts = "results/isoquant/{name}/{name}.transcript_grouped_{token}_counts.tsv".format(name=config["name"], token = GROUP_TOKEN),
            ref_transcript_tpm = "results/isoquant/{name}/{name}.transcript_grouped_{token}_tpm.tsv".format(name=config["name"], token = GROUP_TOKEN),
            ref_gene_counts = "results/isoquant/{name}/{name}.gene_grouped_{token}_counts.tsv".format(name=config["name"], token = GROUP_TOKEN),
            ref_gene_tpm = "results/isoquant/{name}/{name}.gene_grouped_{token}_tpm.tsv".format(name=config["name"], token = GROUP_TOKEN),
            models = "results/isoquant/{name}/{name}.transcript_models.gtf".format(name=config["name"]),
            read2transcripts = "results/isoquant/{name}/{name}.transcript_model_reads.tsv.gz".format(name=config["name"]),
            genedb = GENEDB,
            read_assignments = READ_ASSIGNMENTS
    params: ref = REF,
            basecode_flags = " ".join(flag for flag, on in [
                ("--basecode_correct",            config["basecode_correct"]),
                ("--basecode_keep_nonunique",     config["basecode_keep_nonunique"]),
                ("--basecode_no_context_resolve", config["basecode_no_context_resolve"]),
                ("--basecode_end_resolve",        config["basecode_end_resolve"]),
            ] if on)
    log: "results/isoquant/logs/{name}.run_isoquant.log".format(name=config["name"])
    shell: "python IsoQuant/isoquant.py -p {config[name]} --reference {params.ref} --genedb {input.gff3} --complete_genedb --bam {input.bam} --read_group {config[read_group]} --data_type {config[data_type]} --basecode --basecode_max_gap {config[basecode_max_gap]} {params.basecode_flags} --count_exons --matching_strategy {config[matching_strategy]} --transcript_quantification {config[transcript_quantification]} --gene_quantification {config[gene_quantification]} --model_construction_strategy {config[model_construction_strategy]} --polya_requirement {config[polya_requirement]} --check_canonical --bam_tags TC,IC,FC,NR,ER,IR,AD -o results/isoquant/ --sqanti_output --large_output read_assignments corrected_bed read2transcripts > {log} 2>&1"

rule make_mudata:
    input: disc_transcript_counts = "results/isoquant/{name}/{name}.discovered_transcript_grouped_{token}_counts.linear.tsv".format(name=config["name"], token = GROUP_TOKEN),
           disc_transcript_tpm = "results/isoquant/{name}/{name}.discovered_transcript_grouped_{token}_tpm.tsv".format(name=config["name"], token = GROUP_TOKEN),
           disc_gene_counts = "results/isoquant/{name}/{name}.discovered_gene_grouped_{token}_counts.tsv".format(name=config["name"], token = GROUP_TOKEN),
           disc_gene_tpm = "results/isoquant/{name}/{name}.discovered_gene_grouped_{token}_tpm.tsv".format(name=config["name"], token = GROUP_TOKEN),
           ref_transcript_counts = "results/isoquant/{name}/{name}.transcript_grouped_{token}_counts.tsv".format(name=config["name"], token = GROUP_TOKEN),
           ref_transcript_tpm = "results/isoquant/{name}/{name}.transcript_grouped_{token}_tpm.tsv".format(name=config["name"], token = GROUP_TOKEN),
           ref_gene_counts = "results/isoquant/{name}/{name}.gene_grouped_{token}_counts.tsv".format(name=config["name"], token = GROUP_TOKEN),
           ref_gene_tpm = "results/isoquant/{name}/{name}.gene_grouped_{token}_tpm.tsv".format(name=config["name"], token = GROUP_TOKEN),
           models = "results/isoquant/{name}/{name}.transcript_models.gtf".format(name=config["name"]),
           genedb = GENEDB
    output: "results/isoquant/mudata/{name}_counts.h5mu".format(name=config["name"])
    log: "results/isoquant/logs/{name}.make_mudata.log".format(name=config["name"])
    shell: "python workflow/scripts/make_mudata.py --disc-transcript-counts {input.disc_transcript_counts} --disc-transcript-tpm {input.disc_transcript_tpm} --disc-gene-counts {input.disc_gene_counts} --disc-gene-tpm {input.disc_gene_tpm} --ref-transcript-counts {input.ref_transcript_counts} --ref-transcript-tpm {input.ref_transcript_tpm} --ref-gene-counts {input.ref_gene_counts} --ref-gene-tpm {input.ref_gene_tpm} --models {input.models} --genedb {input.genedb} --output {output} > {log} 2>&1"

rule make_variant_support:
    input: models = "results/isoquant/{name}/{name}.transcript_models.gtf".format(name=config["name"]),
           read2transcripts = "results/isoquant/{name}/{name}.transcript_model_reads.tsv.gz".format(name=config["name"]),
           bam =  "results/isoquant/{name}.isoquant.bam".format(name=config["name"]),
    output: per_variant = "results/isoquant/{name}/{name}.variant_support.per_variant.tsv".format(name=config["name"]),
            summary = "results/isoquant/{name}/{name}.variant_support.summary.txt".format(name=config["name"])
    log: "results/isoquant/logs/{name}.make_variant_support.log".format(name=config["name"])
    shell: "python workflow/scripts/make_variant_support.py --models {input.models} --read2transcripts {input.read2transcripts} --bam {input.bam} > {log} 2>&1"

rule annotate_bam:
    input: bam =  "results/isoquant/{name}.isoquant.bam".format(name=config["name"]),
           bai =  "results/isoquant/{name}.isoquant.bam.bai".format(name=config["name"]),
           read_assignments = READ_ASSIGNMENTS,
           read2transcripts = "results/isoquant/{name}/{name}.transcript_model_reads.tsv.gz".format(name=config["name"])
    output: ANNOTATED_BAM
    log: "results/isoquant/logs/{name}.annotate_bam.log".format(name=config["name"])
    shell: "python workflow/scripts/annotate_bam.py --input {input.bam} --assignments {input.read_assignments} --read2transcripts {input.read2transcripts} --output {output} --duplicate-mode {config[duplicate_mode]} --tag-unassigned > {log} 2>&1"

rule index_annotated_bam:
    input: ANNOTATED_BAM
    output: ANNOTATED_BAM + ".bai"
    log: "results/isoquant/logs/{name}.index_annotated_bam.log".format(name=config["name"])
    shell: "samtools index {input} > {log} 2>&1"