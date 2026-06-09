import argparse
import pysam
import array

def main():
    parser = argparse.ArgumentParser(description='Modify BAM file for IsoQuant', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='Input', type=str, help='Input .bam file')
    parser.add_argument('-o', '--output' ,metavar='Output', type=str, help='Output .bam file')

    args = parser.parse_args()

    bam_infile = args.input
    bam_outfile = args.output
    bam_in = pysam.AlignmentFile(bam_infile, 'rb')
    bam_out = pysam.AlignmentFile(bam_outfile, 'wb', template=bam_in)
    for read in bam_in.fetch(until_eof=True):
        if read.get_tag('TC') > 0 and read.get_tag('FC') > 0:
            read.set_tag('CP', True)
        else:
            read.set_tag('CP', False)
        if read.get_tag('TC') > 0:
            q_q = read.query_qualities
            q_s = read.query_sequence
            cigarstring = read.cigarstring
            if read.is_reverse:
                q_new = array.array('B', 24*[40])
                q_new.extend(q_q)
                s_new = 24*'T'
                s_new = s_new + q_s
                cigarstring_new = '24S' + cigarstring
            else:
                q_new = q_q
                q_new.extend(array.array('B', 24*[40]))
                s_new = q_s
                s_new = s_new + 24*'A'
                cigarstring_new = cigarstring + '24S'
            read.query_sequence = s_new
            read.query_qualities = q_new
            read.cigarstring = cigarstring_new
            bam_out.write(read)
        else:
            bam_out.write(read)


if __name__ == "__main__":
    main()
