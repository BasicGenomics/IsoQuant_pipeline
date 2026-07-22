#!/usr/bin/env python
import argparse
import gzip
import statistics
import pysam

version = "1.0"

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def get_tag_safe(read, tag, default=None):
    try:
        return read.get_tag(tag)
    except KeyError:
        return default

def parse_models(gtf):
    import re
    models = {}
    exon_len = {}
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    gid_re = re.compile(r'gene_id "([^"]+)"')
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            feat = c[2]
            if feat == "transcript":
                m = tid_re.search(c[8])
                if not m:
                    continue
                tid = m.group(1)
                g = gid_re.search(c[8])
                models[tid] = {
                    "chrom": c[0], "strand": c[6], "source": c[1],
                    "gene_id": g.group(1) if g else "",
                    "span": int(c[4]) - int(c[3]) + 1,
                    "is_novel": (c[1] == "IsoQuant"),
                    "n_exons": 0, "model_len": 0,
                }
            elif feat == "exon":
                m = tid_re.search(c[8])
                if not m:
                    continue
                tid = m.group(1)
                exon_len.setdefault(tid, [0, 0])
                exon_len[tid][0] += 1
                exon_len[tid][1] += int(c[4]) - int(c[3]) + 1
    for tid, (n, ln) in exon_len.items():
        if tid in models:
            models[tid]["n_exons"] = n
            models[tid]["model_len"] = ln
    return models

def parse_r2t(path):
    model_reads = {}
    n_star = 0
    with open_maybe_gz(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("read_id\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rid, tid = parts[0], parts[1]
            if tid == "*":
                n_star += 1
                continue
            model_reads.setdefault(tid, []).append(rid)
    return model_reads, n_star

def scan_bam(bam_path, wanted_reads):
    info = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam.fetch(until_eof=True):
        rid = read.query_name
        if rid not in wanted_reads:
            continue
        tc = get_tag_safe(read, "TC", 0) or 0
        fc = get_tag_safe(read, "FC", 0) or 0
        ic = get_tag_safe(read, "IC", 0) or 0
        has_gap = any(op == 2 for op, _ in (read.cigartuples or []))
        cat = ("T" if tc > 0 else "") + ("F" if fc > 0 else "") + ("I" if ic > 0 else "")
        info[rid] = {
            "complete": tc > 0 and fc > 0,      # CP: both ends (full-length)
            "cat": cat or "none",
            "span": read.reference_length or 0,
            "aligned": read.query_alignment_length or 0,
            "has_gap": has_gap,
            "ad": get_tag_safe(read, "AD", ""),
        }
    bam.close()
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", required=True)
    ap.add_argument("--read2transcripts", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--novel-only", action="store_true", help="Report only novel (de-novo, source=IsoQuant) models.")
    args = ap.parse_args()

    r2t = args.read2transcripts
    for suf in (".transcript_model_reads.tsv.gz", ".transcript_model_reads.tsv"):
        if r2t.endswith(suf):
            prefix = r2t[:-len(suf)] + ".variant_support"
            break
    else:
        prefix = r2t + ".variant_support"

    print("Parsing models GTF ...")
    models = parse_models(args.models)
    n_novel = sum(1 for m in models.values() if m["is_novel"])
    print(f"  {len(models):,} models ({n_novel:,} novel / {len(models)-n_novel:,} known)")

    print("Parsing read2transcripts ...")
    model_reads, n_star = parse_r2t(args.read2transcripts)
    wanted = set()
    for tid, reads in model_reads.items():
        if args.novel_only and not models.get(tid, {}).get("is_novel", False):
            continue
        wanted.update(reads)
    print(f"  {len(model_reads):,} models have reads; {len(wanted):,} supporting reads to look up; "
          f"{n_star:,} reads assigned to no model")

    print("Scanning BAM for supporting reads ...")
    info = scan_bam(args.bam, wanted)
    print(f"  found tag/size info for {len(info):,} / {len(wanted):,} reads")

    rows = []
    for tid, reads in model_reads.items():
        m = models.get(tid)
        if m is None:
            continue
        if args.novel_only and not m["is_novel"]:
            continue
        recs = [info[r] for r in reads if r in info]
        n = len(recs)
        if n == 0:
            continue
        spans = [r["span"] for r in recs]
        aligned = [r["aligned"] for r in recs]
        n_complete = sum(1 for r in recs if r["complete"])
        n_gap = sum(1 for r in recs if r["has_gap"])
        cats = {"TFI": 0, "TF": 0, "TI": 0, "FI": 0, "T": 0, "F": 0, "I": 0, "none": 0}
        for r in recs:
            cats[r["cat"]] = cats.get(r["cat"], 0) + 1
        rows.append({
            "transcript_id": tid, "gene_id": m["gene_id"], "chrom": m["chrom"],
            "is_novel": int(m["is_novel"]), "n_exons": m["n_exons"], "model_len": m["model_len"],
            "n_reads": n,
            "n_full_length": n_complete,
            "frac_full_length": round(n_complete / n, 4),
            "n_TF_I": cats["TFI"],
            "n_TF": cats["TF"],
            "n_T_I": cats["TI"],
            "n_F_I": cats["FI"],
            "n_T": cats["T"],
            "n_F": cats["F"],
            "n_I": cats["I"],
            "n_none": cats["none"],
            "mean_span": round(statistics.mean(spans), 1),
            "median_span": int(statistics.median(spans)),
            "mean_aligned": round(statistics.mean(aligned), 1),
            "n_gap": n_gap,
            "frac_gap": round(n_gap / n, 4),
        })

    rows.sort(key=lambda r: r["n_reads"], reverse=True)
    cols = ["transcript_id", "gene_id", "chrom", "is_novel", "n_exons", "model_len",
            "n_reads", "n_full_length", "frac_full_length",
            "n_TF_I", "n_TF", "n_T_I", "n_F_I", "n_T", "n_F", "n_I", "n_none",
            "mean_span", "median_span", "mean_aligned", "n_gap", "frac_gap"]
    out_tsv = prefix + ".per_variant.tsv"
    with open(out_tsv, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

    novel = [r for r in rows if r["is_novel"]]
    known = [r for r in rows if not r["is_novel"]]

    def blk(name, rs):
        if not rs:
            return f"{name}: none\n"
        supp = [r["n_reads"] for r in rs]
        fl = [r["frac_full_length"] for r in rs]
        sz = [r["mean_span"] for r in rs]
        return (f"{name}: {len(rs):,} models\n"
                f"  reads/model:      mean {statistics.mean(supp):.1f}  median {int(statistics.median(supp))}  max {max(supp)}\n"
                f"  full-length frac: mean {statistics.mean(fl):.1%}  median {statistics.median(fl):.1%}\n"
                f"  mean molecule span: median-of-models {int(statistics.median(sz))} bp\n"
                f"  models with >=2 reads: {sum(1 for r in rs if r['n_reads']>=2):,}\n")

    summary = ("=" * 60 + "\nPER-VARIANT SUPPORT SUMMARY\n" + "=" * 60 + "\n"
               + blk("NOVEL (de-novo) models", novel) + "\n"
               + blk("KNOWN models", known)
               + f"\nReads assigned to no model (*): {n_star:,}\n"
               + f"\nWrote: {out_tsv}\n")
    print("\n" + summary)
    with open(prefix + ".summary.txt", "w") as o:
        o.write(summary)


if __name__ == "__main__":
    main()
