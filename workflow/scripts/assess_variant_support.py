#!/usr/bin/env python
import argparse
import gzip
import statistics
import re
import pysam

version = "1.1"

REF_TYPES = ["unique", "unique_minor_difference", "ambiguous",
             "inconsistent_non_intronic", "inconsistent", "inconsistent_ambiguous"]
END_CATS = ["T_F_I", "T_F", "T_I", "F_I", "T", "F", "I", "none"]
CAT_LABEL = {"TFI": "T_F_I", "TF": "T_F", "TI": "T_I", "FI": "F_I",
             "T": "T", "F": "F", "I": "I", "none": "none"}

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def get_tag_safe(read, tag, default=None):
    try:
        return read.get_tag(tag)
    except KeyError:
        return default

def parse_models(gtf):
    models, exon_len = {}, {}
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    gid_re = re.compile(r'gene_id "([^"]+)"')
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            if c[2] == "transcript":
                m = tid_re.search(c[8])
                if not m:
                    continue
                g = gid_re.search(c[8])
                models[m.group(1)] = {"chrom": c[0], "gene_id": g.group(1) if g else "",
                                      "is_novel": (c[1] == "IsoQuant"), "n_exons": 0, "model_len": 0}
            elif c[2] == "exon":
                m = tid_re.search(c[8])
                if not m:
                    continue
                exon_len.setdefault(m.group(1), [0, 0])
                exon_len[m.group(1)][0] += 1
                exon_len[m.group(1)][1] += int(c[4]) - int(c[3]) + 1
    for tid, (n, ln) in exon_len.items():
        if tid in models:
            models[tid]["n_exons"], models[tid]["model_len"] = n, ln
    return models

def parse_r2t(path):
    model_reads, n_star = {}, 0
    with open_maybe_gz(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("read_id\t"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            if p[1] == "*":
                n_star += 1
                continue
            model_reads.setdefault(p[1], []).append(p[0])
    return model_reads, n_star

def parse_read_assignments(path):
    iso_reads, meta, n_unassigned = {}, {}, 0
    with open_maybe_gz(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            rid, chrom, iso, gene, atype = c[0], c[1], c[3], c[4], c[5]
            if iso in (".", "*", ""):
                n_unassigned += 1
                continue
            iso_reads.setdefault(iso, []).append((rid, atype))
            if iso not in meta:
                meta[iso] = (chrom, gene)
    return iso_reads, meta, n_unassigned

def scan_bam(bam_path, wanted):
    info = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam.fetch(until_eof=True):
        rid = read.query_name
        if rid not in wanted:
            continue
        tc = get_tag_safe(read, "TC", 0) or 0
        fc = get_tag_safe(read, "FC", 0) or 0
        ic = get_tag_safe(read, "IC", 0) or 0
        cat = ("T" if tc > 0 else "") + ("F" if fc > 0 else "") + ("I" if ic > 0 else "")
        info[rid] = {
            "complete": tc > 0 and fc > 0,
            "cat": cat or "none",
            "aligned": read.query_alignment_length or 0,
            "has_gap": any(op == 2 for op, _ in (read.cigartuples or [])),
        }
    bam.close()
    return info

def base_row(fid, gene, chrom, recs):
    n = len(recs)
    aligned = [r["aligned"] for r in recs]
    n_complete = sum(1 for r in recs if r["complete"])
    n_gap = sum(1 for r in recs if r["has_gap"])
    cats = {k: 0 for k in CAT_LABEL}
    for r in recs:
        cats[r["cat"]] = cats.get(r["cat"], 0) + 1
    row = {"feature_id": fid, "gene_id": gene, "chrom": chrom, "n_reads": n,
           "n_full_length": n_complete, "frac_full_length": round(n_complete / n, 4)}
    for key, lab in CAT_LABEL.items():
        row["n_" + lab] = cats[key]
    row["mean_aligned"] = round(statistics.mean(aligned), 1)
    row["median_aligned"] = int(statistics.median(aligned))
    row["n_gap"] = n_gap
    row["frac_gap"] = round(n_gap / n, 4)
    return row

def write_tsv(path, rows, cols):
    rows.sort(key=lambda r: r["n_reads"], reverse=True)
    with open(path, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

def summary_block(title, rows, novel_key=None):
    if not rows:
        return f"{title}: none\n"
    supp = [r["n_reads"] for r in rows]
    fl = [r["frac_full_length"] for r in rows]
    out = [f"{title}: {len(rows):,} features",
           f"  reads/feature:    mean {statistics.mean(supp):.1f}  median {int(statistics.median(supp))}  max {max(supp)}",
           f"  full-length frac: mean {statistics.mean(fl):.1%}  median {statistics.median(fl):.1%}",
           f"  features >=2 reads: {sum(1 for r in rows if r['n_reads'] >= 2):,}"]
    if novel_key:
        nv = [r for r in rows if r[novel_key]]
        out.append(f"  novel: {len(nv):,} / known: {len(rows)-len(nv):,}")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["discovered", "reference", "both"], default="both")
    ap.add_argument("--bam", required=True, help="Stitched .bam file")
    ap.add_argument("--prefix", required=True, help="Output path prefix")
    ap.add_argument("--read2transcripts", help="Transcript model reads .tsv.gz file")
    ap.add_argument("--models", help="Transcript models .gtf file")
    ap.add_argument("--read-assignments", help="Read assignments .tsv.gz file")
    args = ap.parse_args()

    do_disc = args.mode in ("discovered", "both")
    do_ref = args.mode in ("reference", "both")
    if do_disc and not (args.read2transcripts and args.models):
        ap.error("--mode discovered/both requires --read2transcripts and --models")
    if do_ref and not args.read_assignments:
        ap.error("--mode reference/both requires --read-assignments")

    wanted = set()
    disc_model_reads = disc_models = None
    ref_iso_reads = ref_meta = None
    n_star = n_unassigned = 0
    if do_disc:
        print("Parsing models GTF + read2transcripts ...")
        disc_models = parse_models(args.models)
        disc_model_reads, n_star = parse_r2t(args.read2transcripts)
        for reads in disc_model_reads.values():
            wanted.update(reads)
    if do_ref:
        print("Parsing read_assignments ...")
        ref_iso_reads, ref_meta, n_unassigned = parse_read_assignments(args.read_assignments)
        for reads in ref_iso_reads.values():
            wanted.update(r for r, _ in reads)

    print(f"Scanning BAM for {len(wanted):,} supporting reads ...")
    info = scan_bam(args.bam, wanted)
    print(f"  found tag/size info for {len(info):,} reads")

    summary_parts = ["=" * 60, "VARIANT / TRANSCRIPT SUPPORT SUMMARY", "=" * 60, ""]

    if do_disc:
        rows = []
        for mid, reads in disc_model_reads.items():
            m = disc_models.get(mid)
            if m is None:
                continue
            recs = [info[r] for r in reads if r in info]
            if not recs:
                continue
            row = base_row(mid, m["gene_id"], m["chrom"], recs)
            row["is_novel"] = int(m["is_novel"])
            row["n_exons"] = m["n_exons"]
            row["model_len"] = m["model_len"]
            rows.append(row)
        cols = (["feature_id", "gene_id", "chrom", "is_novel", "n_exons", "model_len",
                 "n_reads", "n_full_length", "frac_full_length"]
                + ["n_" + c for c in END_CATS]
                + ["mean_aligned", "median_aligned", "n_gap", "frac_gap"])
        out = args.prefix + ".discovered_variant_support.per_variant.tsv"
        write_tsv(out, rows, cols)
        summary_parts.append(summary_block("DISCOVERED models", rows, novel_key="is_novel"))
        summary_parts.append(f"  reads assigned to no model (*): {n_star:,}\n")
        print(f"Wrote {out}")

    if do_ref:
        rows = []
        for iso, reads in ref_iso_reads.items():
            recs = [info[r] for r, _ in reads if r in info]
            if not recs:
                continue
            chrom, gene = ref_meta.get(iso, (".", "."))
            row = base_row(iso, gene, chrom, recs)
            tcount = {t: 0 for t in REF_TYPES}
            tcount["other"] = 0
            for r, atype in reads:
                if r not in info:
                    continue
                if atype in tcount:
                    tcount[atype] += 1
                else:
                    tcount["other"] += 1
            for t in REF_TYPES:
                row["n_" + t] = tcount[t]
            row["n_other_type"] = tcount["other"]
            rows.append(row)
        cols = (["feature_id", "gene_id", "chrom", "n_reads", "n_full_length", "frac_full_length"]
                + ["n_" + t for t in REF_TYPES] + ["n_other_type"]
                + ["n_" + c for c in END_CATS]
                + ["mean_aligned", "median_aligned", "n_gap", "frac_gap"])
        out = args.prefix + ".reference_variant_support.per_variant.tsv"
        write_tsv(out, rows, cols)
        summary_parts.append(summary_block("REFERENCE transcripts", rows))
        summary_parts.append(f"  reads with no reference isoform: {n_unassigned:,}\n")
        print(f"Wrote {out}")

    summary_path = args.prefix + ".variant_support.summary.txt"
    with open(summary_path, "w") as o:
        o.write("\n".join(summary_parts) + "\n")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
