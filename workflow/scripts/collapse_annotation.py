#!/usr/bin/env python3
import argparse, re, sys

version = "1.0"

TX_TYPES = {"transcript", "mRNA"}

def attr(s, k):
    m = re.search(k + r'=([^;]+)', s)
    return m.group(1) if m else ""

def chain(blocks):
    b = sorted(blocks)
    return tuple((b[i][1] + 1, b[i + 1][0] - 1) for i in range(len(b) - 1))

def tag_tier(tags):
    toks = tags.split(",")
    if "MANE_Select" in toks:
        return 0.0
    if "MANE_Plus_Clinical" in toks:
        return 1.0
    if "Ensembl_canonical" in toks:
        return 2.0

    def appris_rank(prefix, base):
        nums, present = [], False
        for t in toks:
            if t.startswith(prefix):
                present = True
                m = re.search(r'_(\d+)$', t)
                if m:
                    nums.append(int(m.group(1)))
        return base + (min(nums) if nums else 9) / 10.0 if present else None

    r = appris_rank("appris_principal", 3.0)
    if r is not None:
        return r
    r = appris_rank("appris_alternative", 4.0)
    if r is not None:
        return r
    if "GENCODE_Primary" in toks:
        return 5.0
    if "CCDS" in toks:
        return 6.0
    if "basic" in toks:
        return 7.0
    return 8.0


def main():
    ap = argparse.ArgumentParser(description="Collapse same-body transcripts; keep most-supported representative")
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--group-by", choices=["chain", "cds", "chain_cds"], default="chain")
    ap.add_argument("--representative", choices=["canonical", "longest"], default="canonical",
                    help="per group keep: 'canonical' = best tag, tie-break longest span; "
                         "'longest' = longest span (longest UTRs) so reads with long UTRs still fit, "
                         "tie-break best tag")
    a = ap.parse_args()

    gene_line = None
    tx = {}
    order = []
    cur = None
    ng = tin = tkept = 0
    fout = open(a.output, "w")

    def key_for(t):
        if a.group_by == "cds":
            return tuple(sorted(t["cds"])) if t["cds"] else ("noCDS",) + chain(t["exons"])
        if a.group_by == "chain_cds":
            return (chain(t["exons"]), tuple(sorted(t["cds"])))
        return chain(t["exons"])

    def structureless(t):
        if a.group_by in ("cds", "chain_cds"):
            return (not t["cds"]) and len(t["exons"]) <= 1
        return len(t["exons"]) <= 1

    def span_iv(t):
        ex = t["exons"]
        return (min(e[0] for e in ex), max(e[1] for e in ex)) if ex else (t["start"], t["end"])

    def flush():
        nonlocal ng, tin, tkept
        if gene_line is None:
            return
        ng += 1
        groups = {}
        loose = []
        for tid in order:
            if structureless(tx[tid]):
                loose.append(tid)
            else:
                groups.setdefault(key_for(tx[tid]), []).append(tid)
        clusters = []
        cur = []
        cmax = None
        for tid in sorted(loose, key=lambda x: span_iv(tx[x])):
            s, e = span_iv(tx[tid])
            if cmax is not None and s <= cmax:
                cur.append(tid); cmax = max(cmax, e)
            else:
                if cur:
                    clusters.append(cur)
                cur = [tid]; cmax = e
        if cur:
            clusters.append(cur)
        keep = set()
        for tids in list(groups.values()) + clusters:
            if a.representative == "longest":
                best = min(tids, key=lambda x: (-tx[x]["span"], tag_tier(tx[x]["tags"])))
            else:
                best = min(tids, key=lambda x: (tag_tier(tx[x]["tags"]), -tx[x]["span"]))
            keep.add(best)
        fout.write(gene_line)
        for tid in order:
            tin += 1
            if tid in keep:
                tkept += 1
                fout.write(tx[tid]["tline"])
                fout.writelines(tx[tid]["children"])
    with open(a.input) as fin:
        for line in fin:
            if line.startswith("#"):
                fout.write(line); continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            typ = c[2]
            if typ == "gene":
                flush(); gene_line = line; tx = {}; order = []; cur = None
            elif typ in TX_TYPES:
                tid = attr(c[8], "transcript_id")
                tx[tid] = {"tline": line, "children": [], "exons": [], "cds": [],
                           "tags": attr(c[8], "tag"), "span": int(c[4]) - int(c[3]),
                           "start": int(c[3]), "end": int(c[4])}
                order.append(tid); cur = tid
            else:
                tid = attr(c[8], "transcript_id") or cur
                if tid in tx:
                    tx[tid]["children"].append(line)
                    if typ == "exon":
                        tx[tid]["exons"].append((int(c[3]), int(c[4])))
                    elif typ == "CDS":
                        tx[tid]["cds"].append((int(c[3]), int(c[4])))
        flush()
    fout.close()
    sys.stderr.write(f"group-by={a.group_by}  genes={ng}  tx_in={tin}  kept={tkept}  "
                     f"dropped={tin-tkept} ({100*(tin-tkept)/max(tin,1):.1f}%)\n")


if __name__ == "__main__":
    main()
