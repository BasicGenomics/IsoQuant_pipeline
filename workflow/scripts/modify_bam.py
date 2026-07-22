import argparse
import pysam
import array

version = "1.1"

def get_tag_safe(read, tag, default=None):
    try:
        return read.get_tag(tag)
    except KeyError:
        return default

def has_deletion(cigartuples):
    return any(op == 2 for op, _ in (cigartuples or []))

def total_deletion(cigartuples):
    return sum(length for op, length in (cigartuples or []) if op == 2)

def main():
    parser = argparse.ArgumentParser(description='Modify BAM file for IsoQuant', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='Input', type=str, help='Input .bam file')
    parser.add_argument('-o', '--output', metavar='Output', type=str, help='Output .bam file')
    parser.add_argument('--full-length-only', action='store_true', help='Emit only full-length molecules (complete: TC>0 & FC>0)')
    parser.add_argument('--tracking-file', metavar='TSV', type=str, default=None, help='Optional .tsv file listing every molecule and what was adapted')

    args = parser.parse_args()

    bam_in = pysam.AlignmentFile(args.input, 'rb')
    bam_out = pysam.AlignmentFile(args.output, 'wb', template=bam_in)

    track = None
    if args.tracking_file:
        track = open(args.tracking_file, 'w')
        track.write("read_id\thas_gap\tdel_len\tpolya_added\tcomplete\tAD\n")

    n_total = n_written = n_dropped_nonfl = 0

    for read in bam_in.fetch(until_eof=True):
        n_total += 1
        tc = get_tag_safe(read, 'TC', 0) or 0
        fc = get_tag_safe(read, 'FC', 0) or 0
        complete = tc > 0 and fc > 0
        read.set_tag('CP', complete)
        has_gap = has_deletion(read.cigartuples)
        del_len = total_deletion(read.cigartuples)
        polya_added = tc > 0
        if polya_added:
            q_q = read.query_qualities
            q_s = read.query_sequence
            cigarstring = read.cigarstring
            if read.is_reverse:
                q_new = array.array('B', 24 * [40])
                q_new.extend(q_q)
                s_new = 24 * 'T' + q_s
                cigarstring_new = '24S' + cigarstring
            else:
                q_new = q_q
                q_new.extend(array.array('B', 24 * [40]))
                s_new = q_s + 24 * 'A'
                cigarstring_new = cigarstring + '24S'
            read.query_sequence = s_new
            read.query_qualities = q_new
            read.cigarstring = cigarstring_new
        if has_gap and polya_added:
            ad = 'gap+polya'
        elif has_gap:
            ad = 'gap'
        elif polya_added:
            ad = 'polya'
        else:
            ad = 'none'
        read.set_tag('AD', ad, value_type='Z')
        if track is not None:
            track.write(f"{read.query_name}\t{int(has_gap)}\t{del_len}\t{int(polya_added)}\t"
                        f"{int(complete)}\t{ad}\n")
        if args.full_length_only and not complete:
            n_dropped_nonfl += 1
            continue
        bam_out.write(read)
        n_written += 1

    bam_in.close()
    bam_out.close()
    if track is not None:
        track.close()
    print(f"Molecules: {n_total:,} read, {n_written:,} written"
          + (f", {n_dropped_nonfl:,} dropped (not full-length)" if args.full_length_only else ""))


if __name__ == "__main__":
    main()
