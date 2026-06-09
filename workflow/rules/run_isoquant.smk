from snakemake.utils import min_version
import time as _bc_time
import sys
import os
min_version("6.7")

configfile: "config/config.yaml"

if 'resource_dir' not in config:
    config['resource_dir'] = '/usr/local/app/resources/'
else:
    if config['resource_dir'] == '':
        config['resource_dir'] = '/usr/local/app/resources/'

if "basecode_max_gap" not in config:
    config["basecode_max_gap"] = 550

if "data_type" not in config:
    config["data_type"] = "pacbio_ccs"

if "stranded" not in config:
    config["stranded"] = "none"

if "model_construction_strategy" not in config:
    config["model_construction_strategy"] = "default_pacbio"

if "matching_strategy" not in config:
    config["matching_strategy"] = "precise"

if "transcript_quantification" not in config:
    config["transcript_quantification"] = "unique_only"

if "gene_quantification" not in config:
    config["gene_quantification"] = "unique_splicing_consistent"

if "polya_requirement" not in config:
    config["polya_requirement"] = "auto"

if "read_group" not in config:
    config["read_group"] = "tag:SM"

GROUP_TOKEN = config["read_group"].replace(":", "_")

REF = "{resource_dir}/genome_references/{reference}/reference.fa".format(resource_dir = config["resource_dir"], reference = config["reference"])
GFF = "{resource_dir}/genome_references/{reference}/geneannotations".format(resource_dir=config["resource_dir"], reference = config["reference"])


rule all:
    input: "results/isoquant/mudata/{name}_counts.h5mu".format(name=config["name"])

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


onstart:
    config['_bc_start'] = _bc_time.time()
    print("▶ BaseCode IsoQuant Pipeline starting", file=sys.stderr, flush=True)
    print(f"  Run ID: {config['name']}  (BaseCode mode, max_gap={config['basecode_max_gap']})", file=sys.stderr, flush=True)

onsuccess:
    elapsed = int(_bc_time.time() - config.get('_bc_start', _bc_time.time()))
    h, rem = divmod(elapsed, 3600)
    m, s = divmod(rem, 60)
    elapsed_str = f"{h}h {m}m {s}s" if h else (f"{m}m {s}s" if m else f"{s}s")
    print(f"✓ BaseCode IsoQuant Pipeline completed in {elapsed_str}", file=sys.stderr, flush=True)

onerror:
    import re as _bc_re
    _e = sys.stderr
    print("", file=_e, flush=True)
    print("✗ BaseCode IsoQuant Pipeline failed", file=_e, flush=True)
    print(f"  Snakemake log: {log}", file=_e, flush=True)
    failed_rules, rule_logs = [], []
    try:
        with open(log) as _fh:
            _txt = _fh.read()
    except Exception:
        _txt = ""
    failed_rules = _bc_re.findall(r"Error in rule (\S+?):", _txt)
    for m in _bc_re.finditer(r"^\s*log:\s*(\S+)", _txt, _bc_re.MULTILINE):
        p = m.group(1).rstrip(",")
        if p not in rule_logs:
            rule_logs.append(p)
    if failed_rules:
        print(f"  Failed rule(s): {', '.join(sorted(set(failed_rules)))}", file=_e, flush=True)
    shown = False
    for p in rule_logs:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            shown = True
            print(f"\n  --- tail of {p} ---", file=_e, flush=True)
            shell("tail -n 40 {p:q} 1>&2")
    if not shown:
        print("\n  --- tail of Snakemake log ---", file=_e, flush=True)
        shell("tail -n 50 {log} 1>&2")
    print("\n  Per-rule logs: results/isoquant/logs/", file=_e, flush=True)