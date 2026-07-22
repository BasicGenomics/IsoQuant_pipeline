#!/usr/bin/env python3

import argparse
import gzip
import sys
import pysam

version = "1.0"

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

def load_assignments(assignments_gz, duplicate_mode="best"):
    assignments = {}
    intern = sys.intern
    priority_lookup = ASSIGNMENT_PRIORITY
    with gzip.open(assignments_gz, "rt") as f:
        header = None
        for line in f:
            if line.lstrip("#").startswith("read_id\t"):
                header = [h.lstrip("#") for h in line.rstrip("\n").split("\t")]
                break
        if header is None:
            raise ValueError("Could not find read_id header line")
        try:
            rid_i = header.index("read_id")
            atype_i = header.index("assignment_type")
            events_i = header.index("assignment_events")
            iso_i = header.index("isoform_id")
            gene_i = header.index("gene_id")
        except ValueError as e:
            raise ValueError(f"Missing required column: {e}")
        rows_seen = 0
        for line in f:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            rid = fields[rid_i]
            atype = fields[atype_i]
            entry = (
                intern(atype),
                intern(fields[events_i]),
                intern(fields[iso_i]),
                intern(fields[gene_i]),
            )
            if duplicate_mode == "all":
                assignments.setdefault(rid, []).append(entry)
            elif duplicate_mode == "first":
                if rid in assignments:
                    continue
                assignments[rid] = entry
            else:
                existing = assignments.get(rid)
                if existing is not None and priority_lookup.get(existing[0], 999) <= priority_lookup.get(atype, 999):
                    continue
                assignments[rid] = entry
            rows_seen += 1
            if rows_seen % 1_000_000 == 0:
                print(f"  Read {rows_seen:,} assignment rows ({len(assignments):,} unique read_ids)")
    return assignments


def load_read2transcripts(path):
    r2t = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#") or line.startswith("read_id\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rid, model = parts[0], parts[1]
            if model == "*":
                continue
            if rid in r2t:
                r2t[rid] += ";" + model
            else:
                r2t[rid] = model
    return r2t


def tag_bam(
    input_bam,
    output_bam,
    assignments_gz,
    duplicate_mode="best",
    tag_unassigned=False,
    read2transcripts=None,
):
    assignments = load_assignments(assignments_gz, duplicate_mode)
    print(f"Loaded assignments for {len(assignments):,} read_ids")
    r2t = {}
    if read2transcripts:
        r2t = load_read2transcripts(read2transcripts)
        print(f"Loaded discovered-model assignments for {len(r2t):,} read_ids")
    total = tagged = missing = 0
    with pysam.AlignmentFile(input_bam, "rb") as bam_in, \
         pysam.AlignmentFile(output_bam, "wb", header=bam_in.header) as bam_out:
        for read in bam_in:
            total += 1
            entry = assignments.get(read.query_name)
            if entry is not None:
                if duplicate_mode == "all":
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
                tagged += 1
            else:
                missing += 1
                if tag_unassigned:
                    read.set_tag("ZA", "no_isoquant_assignment", "Z", replace=True)
            if read2transcripts:
                model = r2t.get(read.query_name)
                if model is not None:
                    read.set_tag("ZM", model, "Z", replace=True)
                elif tag_unassigned:
                    read.set_tag("ZM", "no_model", "Z", replace=True)
            bam_out.write(read)
            if total % 1_000_000 == 0:
                print(f"Processed {total:,} | tagged {tagged:,} | missing {missing:,}")
    print(f"Total BAM records processed: {total:,}")
    print(f"Tagged BAM records: {tagged:,}")
    print(f"BAM records without IsoQuant assignment: {missing:,}")


def main():
    parser = argparse.ArgumentParser(description="Add IsoQuant read assignment tags to BAM file")
    parser.add_argument("-i", "--input", required=True, help="Input .bam file")
    parser.add_argument("-a", "--assignments", required=True, help="Read assignments .tsv.gz file")
    parser.add_argument("-o", "--output", required=True, help="Output .bam file")
    parser.add_argument("--read2transcripts", default=None, help="Transcript model reads .tsv.gz file")
    parser.add_argument("--duplicate-mode", choices=["first", "best", "all"], default="all",
                        help="Resolve multi-isoform rows per read: "
                             "'first' = first row; "
                             "'best' = highest-priority type (1 isoform); "
                             "'all' = every compatible isoform, ZI and ZE ';'-joined and positionally aligned (ZI[i] <-> ZE[i])")
    parser.add_argument("--tag-unassigned", action="store_true", help="Tag reads without assignment")
    args = parser.parse_args()
 
    tag_bam(
        input_bam=args.input,
        output_bam=args.output,
        assignments_gz=args.assignments,
        duplicate_mode=args.duplicate_mode,
        tag_unassigned=args.tag_unassigned,
        read2transcripts=args.read2transcripts,
    )


if __name__ == "__main__":
    main()
