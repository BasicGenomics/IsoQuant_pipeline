#!/usr/bin/env python3

import argparse
import array
import gzip
import os
import re
import shutil
import sys
import pysam

version = "1.2.0"

DEFAULT_DROP_TAGS = ("CV",)

ASSIGNMENT_PRIORITY = {
    "unique": 1,
    "unique_minor_difference": 2,
    "ambiguous": 3,
    "inconsistent_ambiguous": 4,
    "inconsistent_non_intronic": 5,
    "inconsistent": 6,
    "noninformative": 7,
    "intergenic": 8,
}

def _open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


class AssignmentReader:

    def __init__(self, path, duplicate_mode="best"):
        self.f = _open_text(path)
        self.duplicate_mode = duplicate_mode
        self._pushback = None
        self.exhausted = False
        self.rows_seen = 0
        header = None
        for line in self.f:
            if line.lstrip("#").startswith("read_id\t"):
                header = [h.lstrip("#") for h in line.rstrip("\n").split("\t")]
                break
        if header is None:
            raise ValueError("Could not find read_id header line")
        try:
            self.rid_i = header.index("read_id")
            self.chr_i = header.index("chr")
            self.atype_i = header.index("assignment_type")
            self.events_i = header.index("assignment_events")
            self.iso_i = header.index("isoform_id")
            self.gene_i = header.index("gene_id")
        except ValueError as e:
            raise ValueError(f"Missing required column: {e}")

    def _next_row(self):
        if self._pushback is not None:
            row, self._pushback = self._pushback, None
            return row
        for line in self.f:
            if not line or line[0] == "#":
                continue
            return line.rstrip("\n").split("\t")
        self.exhausted = True
        return None

    def next_block(self):
        row = self._next_row()
        if row is None:
            return None, None
        chrom = row[self.chr_i]
        intern = sys.intern
        priority = ASSIGNMENT_PRIORITY
        mode = self.duplicate_mode
        rid_i, atype_i, events_i, iso_i, gene_i = (
            self.rid_i, self.atype_i, self.events_i, self.iso_i, self.gene_i)
        out = {}
        while row is not None:
            if row[self.chr_i] != chrom:
                self._pushback = row
                break
            rid = row[rid_i]
            atype = row[atype_i]
            entry = (intern(atype), intern(row[events_i]),
                     intern(row[iso_i]), intern(row[gene_i]))
            if mode == "all":
                existing = out.get(rid)
                if existing is None:
                    out[rid] = [entry]
                else:
                    existing.append(entry)
            elif mode == "first":
                if rid not in out:
                    out[rid] = entry
            else:
                existing = out.get(rid)
                if existing is None or priority.get(existing[0], 999) > priority.get(atype, 999):
                    out[rid] = entry
            self.rows_seen += 1
            row = self._next_row()
        return chrom, out


class BedBlockReader:

    def __init__(self, path):
        self.f = _open_text(path)
        self.stash = {}
        self._pushback = None
        self.exhausted = False
        self.rows_seen = 0

    def _next_row(self):
        if self._pushback is not None:
            row, self._pushback = self._pushback, None
            return row
        for line in self.f:
            if not line or line[0] == "#" or line.startswith("track") or not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 12:
                continue
            return row
        self.exhausted = True
        return None

    @staticmethod
    def _parse(row):
        cstart = int(row[1])
        sizes = [int(x) for x in row[10].rstrip(",").split(",") if x != ""]
        starts = [int(x) for x in row[11].rstrip(",").split(",") if x != ""]
        return cstart, [(cstart + starts[i], sizes[i]) for i in range(len(sizes))]

    def take(self, chrom):
        out = self.stash.pop(chrom, {})
        if self.exhausted:
            return out
        row = self._next_row()
        while row is not None:
            c = row[0]
            if c == chrom:
                out[row[3]] = self._parse(row)
                self.rows_seen += 1
            elif not out and not self.stash:
                self.stash.setdefault(c, {})[row[3]] = self._parse(row)
                self.rows_seen += 1
            else:
                self._pushback = row
                break
            row = self._next_row()
        return out


class ModelReadsReader:

    MISS_RUN = 4096

    def __init__(self, path):
        self.f = _open_text(path)
        self.carry = []
        self.exhausted = False
        self.rows_seen = 0

    def _next_row(self):
        for line in self.f:
            if not line or line[0] == "#" or line.startswith("read_id\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[1] == "*":
                continue
            return parts[0], parts[1]
        self.exhausted = True
        return None

    @staticmethod
    def _add(out, rid, model):
        existing = out.get(rid)
        out[rid] = model if existing is None else existing + ";" + model

    def take(self, wanted):
        out = {}
        if self.carry:
            keep = []
            for rid, model in self.carry:
                if rid in wanted:
                    self._add(out, rid, model)
                else:
                    keep.append((rid, model))
            self.carry = keep
        if self.exhausted or not wanted:
            return out
        intern = sys.intern
        misses = 0
        while True:
            row = self._next_row()
            if row is None:
                break
            rid, model = row
            if rid in wanted:
                self._add(out, rid, intern(model))
                self.rows_seen += 1
                misses = 0
            else:
                self.carry.append((rid, intern(model)))
                misses += 1
                if misses >= self.MISS_RUN:
                    break
        return out


def cigar_from_blocks(exons):
    cig, prev_end = [], None
    for gstart, size in exons:
        if prev_end is not None and gstart - prev_end > 0:
            cig.append((3, gstart - prev_end))
        cig.append((0, size))
        prev_end = gstart + size
    return cig


class Tagger:

    def __init__(self, duplicate_mode, tag_unassigned, drop_tags, use_bed):
        self.duplicate_mode = duplicate_mode
        self.tag_unassigned = tag_unassigned
        self.drop_tags = drop_tags
        self.use_bed = use_bed
        self.assignments = {}
        self.blocks = {}
        self.r2t = {}
        self.use_r2t = False
        self.total = self.tagged = self.missing = self.imputed = self.dropped = 0

    def process(self, read, out):
        self.total += 1
        for t in self.drop_tags:
            if read.has_tag(t):
                read.set_tag(t, None)
                self.dropped += 1
        if self.use_bed:
            bentry = self.blocks.get(read.query_name)
            if bentry is not None:
                orig_cigar = read.cigarstring or ""
                cstart, exons = bentry
                new_cig = cigar_from_blocks(exons)
                seqlen = sum(sz for _, sz in exons)
                read.query_sequence = "N" * seqlen
                read.reference_start = cstart
                read.cigartuples = new_cig
                read.query_qualities = array.array("B", [30] * seqlen)
                read.set_tag("OC", orig_cigar, "Z", replace=True)
                changed = read.cigarstring != orig_cigar
                read.set_tag("IM", 1 if changed else 0, "i", replace=True)
                if changed:
                    self.imputed += 1
            else:
                read.set_tag("IM", 0, "i", replace=True)
        entry = self.assignments.get(read.query_name)
        if entry is not None:
            if self.duplicate_mode == "all":
                atypes = [e[0] for e in entry]
                genes = [e[3] for e in entry]
                atype = atypes[0] if len(set(atypes)) == 1 else ";".join(atypes)
                gene = genes[0] if len(set(genes)) == 1 else ";".join(dict.fromkeys(genes))
                events = ";".join(e[1] for e in entry)
                iso = ";".join(e[2] for e in entry)
            else:
                atype, events, iso, gene = entry
            read.set_tag("ZA", atype, "Z", replace=True)
            read.set_tag("ZE", events, "Z", replace=True)
            read.set_tag("ZI", iso, "Z", replace=True)
            read.set_tag("ZG", gene, "Z", replace=True)
            self.tagged += 1
        else:
            self.missing += 1
            if self.tag_unassigned:
                read.set_tag("ZA", "no_isoquant_assignment", "Z", replace=True)
        if self.use_r2t:
            model = self.r2t.get(read.query_name)
            if model is not None:
                read.set_tag("ZM", model, "Z", replace=True)
            elif self.tag_unassigned:
                read.set_tag("ZM", "no_model", "Z", replace=True)
        out.write(read)


class PartWriter:
    def __init__(self, path, key, header, threads):
        self.path = path
        self.key = key
        self.header = header
        self.threads = threads
        self.fh = None
        self.n = 0
        self.ordered = True
        self._prev = -1

    def write(self, read):
        if self.fh is None:
            self.fh = pysam.AlignmentFile(self.path, "wb", header=self.header,
                                          threads=self.threads)
        pos = read.reference_start
        if pos < self._prev:
            self.ordered = False
        else:
            self._prev = pos
        self.fh.write(read)
        self.n += 1

    def close(self):
        if self.fh is not None:
            self.fh.close()
        return self.n


def merge_parts(parts, output_bam, part_dir, threads):
    parts = [p for p in sorted(parts, key=lambda x: x.key) if p.n]
    if not parts:
        sys.exit("ERROR: no records written; refusing to produce an empty BAM")
    unsorted = [p for p in parts if not p.ordered]
    if unsorted:
        print(f"Sorting {len(unsorted)} of {len(parts)} chromosome parts "
              f"(reads moved by the imputed-structure rewrite)")
    else:
        print(f"All {len(parts)} chromosome parts came out coordinate-sorted; no sort needed")
    final = []
    for part in parts:
        if part.ordered:
            final.append(part.path)
            continue
        out = part.path[:-4] + ".sorted.bam"
        pysam.sort("-@", str(threads), "-o", out, part.path)
        os.remove(part.path)
        final.append(out)
    print(f"Concatenating {len(final)} parts in BAM header order -> {output_bam}")
    if len(final) > 500:
        fofn = os.path.join(part_dir, "parts.fofn")
        with open(fofn, "w") as fh:
            fh.write("\n".join(final) + "\n")
        pysam.cat("-o", output_bam, "-b", fofn)
    else:
        pysam.cat("-o", output_bam, *final)
    pysam.index("-@", str(threads), output_bam)
    shutil.rmtree(part_dir, ignore_errors=True)

def tag_bam(
    input_bam,
    output_bam,
    assignments_gz,
    duplicate_mode="best",
    tag_unassigned=False,
    read2transcripts=None,
    corrected_bed=None,
    drop_tags=DEFAULT_DROP_TAGS,
    threads=4,
):
    drop_tags = tuple(drop_tags or ())
    if drop_tags:
        print(f"Dropping input tags: {', '.join(drop_tags)}")
    print(f"BGZF threads: {threads}")

    areader = AssignmentReader(assignments_gz, duplicate_mode)
    breader = BedBlockReader(corrected_bed) if corrected_bed else None
    mreader = ModelReadsReader(read2transcripts) if read2transcripts else None

    tagger = Tagger(duplicate_mode, tag_unassigned, drop_tags, bool(corrected_bed))
    tagger.use_r2t = bool(read2transcripts)

    part_dir = output_bam + ".parts"
    shutil.rmtree(part_dir, ignore_errors=True)
    os.makedirs(part_dir, exist_ok=True)
    parts = []
    UNMAPPED_KEY = 10 ** 9

    seen_refs = set()
    with pysam.AlignmentFile(input_bam, "rb", threads=threads) as bam_in:
        ref_index = {r: i for i, r in enumerate(bam_in.references)}

        def new_part(key, label):
            safe = re.sub(r"[^A-Za-z0-9._-]", "_", label)
            part = PartWriter(os.path.join(part_dir, "%010d_%s.bam" % (key, safe)),
                              key, bam_in.header, threads)
            parts.append(part)
            return part

        while True:
            chrom, assignments = areader.next_block()
            if chrom is None:
                break
            tagger.assignments = assignments
            tagger.blocks = breader.take(chrom) if breader else {}
            if mreader:
                wanted = set(assignments)
                if tagger.blocks:
                    wanted.update(tagger.blocks)
                tagger.r2t = mreader.take(wanted)
            if chrom not in ref_index:
                print(f"  WARNING: '{chrom}' is in the assignments file but not in the BAM header "
                      f"- {len(assignments):,} assignments unused (chromosome naming mismatch?)")
                continue
            if chrom in seen_refs:
                sys.exit(f"ERROR: '{chrom}' appears in more than one block of {assignments_gz}. "
                         f"Each block fetches the whole chromosome, so every record on '{chrom}' "
                         f"would be written twice. Group the assignments file by chromosome.")
            seen_refs.add(chrom)
            part = new_part(ref_index[chrom], chrom)
            for read in bam_in.fetch(chrom):
                tagger.process(read, part)
            part.close()
            print(f"  {chrom}: {len(assignments):,} assigned reads | "
                  f"running total {tagger.total:,} records, {tagger.tagged:,} tagged")
            tagger.assignments = tagger.blocks = tagger.r2t = {}

        tagger.assignments = tagger.blocks = tagger.r2t = {}
        n_before = tagger.total
        for ref in bam_in.references:
            if ref in seen_refs:
                continue
            part = new_part(ref_index[ref], ref)
            for read in bam_in.fetch(ref):
                tagger.process(read, part)
            part.close()
        try:
            part = new_part(UNMAPPED_KEY, "unmapped")
            for read in bam_in.fetch("*"):
                tagger.process(read, part)
            part.close()
        except ValueError:
            pass
        if tagger.total > n_before:
            print(f"  unassigned references / unmapped: {tagger.total - n_before:,} records")

    print(f"Total BAM records processed: {tagger.total:,}")
    print(f"Tagged BAM records: {tagger.tagged:,}")
    print(f"BAM records without IsoQuant assignment: {tagger.missing:,}")
    print(f"Assignment rows read: {areader.rows_seen:,}")
    if breader:
        print(f"Reads with CIGAR rewritten to imputed structure (IM=1): {tagger.imputed:,}")
        print(f"corrected_bed rows read: {breader.rows_seen:,}")
    if mreader:
        print(f"Discovered-model rows applied: {mreader.rows_seen:,}")
        if mreader.carry:
            print(f"  NOTE: {len(mreader.carry):,} model rows matched no read on any chromosome")
    if drop_tags and tagger.dropped:
        print(f"Stripped {', '.join(drop_tags)} from {tagger.dropped:,} records")

    merge_parts(parts, output_bam, part_dir, threads)

def main():
    parser = argparse.ArgumentParser(description="Add IsoQuant read assignment tags to BAM file")
    parser.add_argument("-i", "--input", required=True, help="Input .bam file (coordinate-sorted and indexed)")
    parser.add_argument("-a", "--assignments", required=True, help="Read assignments .tsv.gz file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output .bam file; written coordinate-sorted and indexed")
    parser.add_argument("--read2transcripts", default=None, help="Transcript model reads .tsv.gz file")
    parser.add_argument("--duplicate-mode", choices=["first", "best", "all"], default="all",
                        help="Resolve multi-isoform rows per read: "
                             "'first' = first row; "
                             "'best' = highest-priority type (1 isoform); "
                             "'all' = every compatible isoform, ZI and ZE ';'-joined and positionally aligned (ZI[i] <-> ZE[i])")
    parser.add_argument("--tag-unassigned", action="store_true", help="Tag reads without assignment")
    parser.add_argument("--drop-tags", default=",".join(DEFAULT_DROP_TAGS),
                        help="Comma-separated tags to strip from the input records "
                             "(default: %(default)s; pass '' to keep everything)")
    parser.add_argument("--corrected-bed", default=None,
                        help="IsoQuant corrected_reads.bed[.gz]; if given, also rewrite each read's "
                             "CIGAR to its imputed exon/intron structure (SEQ->N, OC=orig cigar, IM=1 if changed)")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="BGZF compression/decompression threads for the input and output BAM "
                             "(default: %(default)s)")
    args = parser.parse_args()

    if not (os.path.exists(args.input + ".bai") or os.path.exists(args.input[:-4] + ".bai")):
        sys.exit(f"ERROR: no .bai next to {args.input}; annotate_bam {version} reads the BAM one "
                 f"chromosome at a time and needs the index. Run 'samtools index' first.")

    tag_bam(
        input_bam=args.input,
        output_bam=args.output,
        assignments_gz=args.assignments,
        duplicate_mode=args.duplicate_mode,
        tag_unassigned=args.tag_unassigned,
        read2transcripts=args.read2transcripts,
        corrected_bed=args.corrected_bed,
        drop_tags=[t.strip() for t in args.drop_tags.split(",") if t.strip()],
        threads=args.threads,
    )


if __name__ == "__main__":
    main()
