#!/usr/bin/env python
import argparse
import gzip
import os
import statistics
import re
import sys
import pysam

version = "1.1.0"

REF_TYPES = ["unique", "unique_minor_difference", "ambiguous",
             "inconsistent_non_intronic", "inconsistent", "inconsistent_ambiguous"]
END_CATS = ["T_F_I", "T_F", "T_I", "F_I", "T", "F", "I", "none"]
CAT_LABEL = {"TFI": "T_F_I", "TF": "T_F", "TI": "T_I", "FI": "F_I",
             "T": "T", "F": "F", "I": "I", "none": "none"}

I_COMPLETE, I_CAT, I_ALIGNED, I_GAP = 0, 1, 2, 3

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
                models[m.group(1)] = {"chrom": sys.intern(c[0]),
                                      "gene_id": g.group(1) if g else "",
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

def assignment_blocks(path, counter):
    with open_maybe_gz(path) as f:
        chrom, iso_reads, meta = None, {}, {}
        for line in f:
            if line.startswith("#") or line.lstrip("#").startswith("read_id\t"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            rid, rchrom, iso, gene, atype = c[0], c[1], c[3], c[4], c[5]
            if rchrom != chrom:
                if chrom is not None:
                    yield chrom, iso_reads, meta
                chrom, iso_reads, meta = rchrom, {}, {}
            if iso in (".", "*", ""):
                counter["unassigned"] += 1
                continue
            iso_reads.setdefault(iso, []).append((rid, sys.intern(atype)))
            if iso not in meta:
                meta[iso] = (sys.intern(rchrom), gene)
        if chrom is not None:
            yield chrom, iso_reads, meta


class ModelReadsByChrom:
    SCAN_LIMIT = 2_000_000

    def __init__(self, path, tid_chrom):
        self.f = open_maybe_gz(path)
        self.tid_chrom = tid_chrom
        self.stash = {}
        self._pushback = None
        self.exhausted = False
        self.n_star = 0
        self.n_unknown_tid = 0

    def _next_row(self):
        if self._pushback is not None:
            row, self._pushback = self._pushback, None
            return row
        for line in self.f:
            if line.startswith("#") or line.startswith("read_id\t"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            if p[1] == "*":
                self.n_star += 1
                continue
            return p[0], p[1]
        self.exhausted = True
        return None

    def take(self, chrom):
        out = self.stash.pop(chrom, {})
        if self.exhausted:
            return out
        started = bool(out)
        scanned = 0
        while True:
            row = self._next_row()
            if row is None:
                break
            rid, tid = row
            c = self.tid_chrom.get(tid)
            if c is None:
                self.n_unknown_tid += 1
                continue
            if c == chrom:
                out.setdefault(tid, []).append(rid)
                started = True
            elif started:
                self._pushback = (rid, tid)
                break
            else:
                self.stash.setdefault(c, {}).setdefault(tid, []).append(rid)
                scanned += 1
                if scanned >= self.SCAN_LIMIT:
                    break
        return out


def scan_bam_chrom(bam, chrom, wanted):
    info = {}
    if not wanted:
        return info
    intern = sys.intern
    for read in bam.fetch(chrom):
        rid = read.query_name
        if rid not in wanted:
            continue
        tc = get_tag_safe(read, "TC", 0) or 0
        fc = get_tag_safe(read, "FC", 0) or 0
        ic = get_tag_safe(read, "IC", 0) or 0
        cat = ("T" if tc > 0 else "") + ("F" if fc > 0 else "") + ("I" if ic > 0 else "")
        info[rid] = (tc > 0 and fc > 0,
                     intern(cat or "none"),
                     read.query_alignment_length or 0,
                     any(op == 2 for op, _ in (read.cigartuples or [])))
    return info

def base_row(fid, gene, chrom, recs):
    n = len(recs)
    aligned = [r[I_ALIGNED] for r in recs]
    n_complete = sum(1 for r in recs if r[I_COMPLETE])
    n_gap = sum(1 for r in recs if r[I_GAP])
    cats = {k: 0 for k in CAT_LABEL}
    for r in recs:
        cats[r[I_CAT]] = cats.get(r[I_CAT], 0) + 1
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
    ap.add_argument("--bam", required=True, help="Stitched .bam file (coordinate-sorted and indexed)")
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
    if not args.read_assignments:
        ap.error("chromosome-by-chromosome processing needs --read-assignments to drive it")
    if not (os.path.exists(args.bam + ".bai") or os.path.exists(args.bam[:-4] + ".bai")):
        sys.exit(f"ERROR: no .bai next to {args.bam}; assess_variant_support {version} reads the BAM "
                 f"one chromosome at a time and needs the index. Run 'samtools index' first.")

    disc_models = mreader = None
    if do_disc:
        print("Parsing models GTF ...")
        disc_models = parse_models(args.models)
        print(f"  {len(disc_models):,} transcript models")
        mreader = ModelReadsByChrom(args.read2transcripts,
                                    {t: m["chrom"] for t, m in disc_models.items()})

    counter = {"unassigned": 0}
    disc_rows, ref_rows = [], []
    n_reads_seen = 0

    bam = pysam.AlignmentFile(args.bam, "rb")
    refs = set(bam.references)
    for chrom, ref_iso_reads, ref_meta in assignment_blocks(args.read_assignments, counter):
        disc_model_reads = mreader.take(chrom) if mreader else {}
        wanted = set()
        if do_ref:
            for reads in ref_iso_reads.values():
                wanted.update(r for r, _ in reads)
        if do_disc:
            for reads in disc_model_reads.values():
                wanted.update(reads)
        if chrom not in refs:
            print(f"  WARNING: '{chrom}' not in the BAM header - {len(wanted):,} reads skipped "
                  f"(chromosome naming mismatch?)")
            continue
        info = scan_bam_chrom(bam, chrom, wanted)
        n_reads_seen += len(info)

        if do_disc:
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
                disc_rows.append(row)

        if do_ref:
            for iso, reads in ref_iso_reads.items():
                recs = [info[r] for r, _ in reads if r in info]
                if not recs:
                    continue
                rchrom, gene = ref_meta.get(iso, (".", "."))
                row = base_row(iso, gene, rchrom, recs)
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
                ref_rows.append(row)

        print(f"  {chrom}: {len(wanted):,} reads | {len(info):,} with tag info | "
              f"running rows disc={len(disc_rows):,} ref={len(ref_rows):,}")
    bam.close()

    if mreader and mreader.stash:
        left = sum(len(v) for v in mreader.stash.values())
        print(f"  NOTE: {left:,} transcript models had reads but no assignment block; not reported")

    print(f"Reads with tag/size info: {n_reads_seen:,}")
    summary_parts = ["=" * 60, "VARIANT / TRANSCRIPT SUPPORT SUMMARY", "=" * 60, ""]

    if do_disc:
        cols = (["feature_id", "gene_id", "chrom", "is_novel", "n_exons", "model_len",
                 "n_reads", "n_full_length", "frac_full_length"]
                + ["n_" + c for c in END_CATS]
                + ["mean_aligned", "median_aligned", "n_gap", "frac_gap"])
        out = args.prefix + ".discovered_variant_support.per_variant.tsv"
        write_tsv(out, disc_rows, cols)
        summary_parts.append(summary_block("DISCOVERED models", disc_rows, novel_key="is_novel"))
        summary_parts.append(f"  reads assigned to no model (*): {mreader.n_star:,}\n")
        print(f"Wrote {out}")

    if do_ref:
        cols = (["feature_id", "gene_id", "chrom", "n_reads", "n_full_length", "frac_full_length"]
                + ["n_" + t for t in REF_TYPES] + ["n_other_type"]
                + ["n_" + c for c in END_CATS]
                + ["mean_aligned", "median_aligned", "n_gap", "frac_gap"])
        out = args.prefix + ".reference_variant_support.per_variant.tsv"
        write_tsv(out, ref_rows, cols)
        summary_parts.append(summary_block("REFERENCE transcripts", ref_rows))
        summary_parts.append(f"  reads with no reference isoform: {counter['unassigned']:,}\n")
        print(f"Wrote {out}")

    summary_path = args.prefix + ".variant_support.summary.txt"
    with open(summary_path, "w") as o:
        o.write("\n".join(summary_parts) + "\n")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
