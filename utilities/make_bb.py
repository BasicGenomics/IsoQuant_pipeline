#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

def run(cmd, shell=False, check=True):
    """Run a command, print it, and optionally check return code."""
    print(f"[cmd] {cmd if isinstance(cmd, str) else ' '.join(cmd)}")
    result = subprocess.run(cmd, shell=shell, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")
    if check and result.returncode != 0:
        sys.exit(result.returncode)
    return result

def need(tool):
    if shutil.which(tool) is None:
        print(f"ERROR: required tool not found in PATH: {tool}", file=sys.stderr)
        sys.exit(127)

def main():
    p = argparse.ArgumentParser(
        description="Generate bigBed (bigGenePred) from GTF/GFF3",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--gtf", required=True, help="GTF or GFF3 annotation file (.gtf/.gff3[.gz])")
    p.add_argument("--fasta", required=True, help="Genome FASTA file")
    p.add_argument("-o", "--output", required=True, help="Output .bb file")
    p.add_argument("--no_download_as", action="store_true",
                   help="Do NOT auto-download bigGenePred.as (you must provide one and pass --as_path)")
    p.add_argument("--as_path", default=None, help="Path to bigGenePred.as (overrides auto-download)")
    args = p.parse_args()

    # Required tools
    need("bedToBigBed")
    need("genePredToBigGenePred")
    need("gtfToGenePred")
    need("gff3ToGenePred")
    need("samtools")
    need("sort")
    need("awk")

    base = os.path.splitext(os.path.basename(args.gtf))[0]
    # Handle .gtf.gz / .gff3.gz
    is_gz = args.gtf.endswith(".gz")
    is_gtf = args.gtf.endswith(".gtf") or args.gtf.endswith(".gtf.gz")
    is_gff3 = args.gtf.endswith(".gff3") or args.gtf.endswith(".gff3.gz")
    if not (is_gtf or is_gff3):
        print("ERROR: --gtf must be .gtf(.gz) or .gff3(.gz)", file=sys.stderr)
        sys.exit(2)

    work = tempfile.mkdtemp(prefix="bbuild_")
    gene_pred = os.path.join(work, f"{base}.genePred")
    biggene_txt = os.path.join(work, f"{base}.bigGenePred.txt")
    biggene_sorted = os.path.join(work, f"{base}.bigGenePred.sorted.txt")
    chrom_sizes = os.path.join(work, f"{base}.chrom.sizes")

    # 1) GTF/GFF3 -> genePred (accept .gz via zcat)
    if is_gtf:
        if is_gz:
            run(f"zcat {args.gtf} | gtfToGenePred -genePredExt /dev/stdin {gene_pred}", shell=True)
        else:
            run(["gtfToGenePred", "-genePredExt", args.gtf, gene_pred])
    else:
        if is_gz:
            run(f"zcat {args.gtf} | gff3ToGenePred /dev/stdin {gene_pred}", shell=True)
        else:
            run(["gff3ToGenePred", args.gtf, gene_pred])

    # 2) genePred -> bigGenePred (tab text)
    run(["genePredToBigGenePred", gene_pred, biggene_txt])

    # 3) Ensure sorted by chrom, chromStart (required by bedToBigBed)
    # bigGenePred text fields start with: chrom, chromStart, chromEnd, ...
    run(f"LC_ALL=C sort -k1,1 -k2,2n {biggene_txt} > {biggene_sorted}", shell=True)

    # 4) chrom.sizes from FASTA index (build if missing)
    fai = f"{args.fasta}.fai"
    if not os.path.exists(fai):
        run(["samtools", "faidx", args.fasta])
    run(f"awk '{{print $1\"\\t\"$2}}' {fai} > {chrom_sizes}", shell=True)

    # 5) Get/choose .as schema (bigGenePred)
    as_path = args.as_path
    if as_path is None and not args.no_download_as:
        as_path = os.path.join(work, "bigGenePred.as")
        # Official schema from UCSC
        run(f"wget -qO {as_path} https://genome.ucsc.edu/goldenPath/help/examples/bigGenePred.as", shell=True)
        if not os.path.exists(as_path) or os.path.getsize(as_path) == 0:
            print("ERROR: failed to download bigGenePred.as. Provide --as_path or try again.", file=sys.stderr)
            sys.exit(3)

    # 6) Make .bb
    # If you provide a .as, we treat it as bed12+8 (canonical bigGenePred)
    if as_path:
        run([
            "bedToBigBed",
            "-type=bed12+8",
            "-tab",
            f"-as={as_path}",
            biggene_sorted,
            chrom_sizes,
            args.output
        ])
    else:
        # Technically possible to make bigBed without .as, but not recommended
        run([
            "bedToBigBed",
            biggene_sorted,
            chrom_sizes,
            args.output
        ])

    print(f" bigBed path: {args.output}")

if __name__ == "__main__":
    main()
