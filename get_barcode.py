import argparse
import os


def get_counts(sam_file, fp_len, fp_start):
    import pysam
    import pandas as pd
    from collections import Counter

    sam = pysam.Samfile(sam_file)
    bcs = Counter([x.seq[(fp_start - x.reference_start - 1):][:fp_len] for x in sam if
                   x.reference_start <= fp_start and x.reference_end >= fp_start + fp_len])
    for key in list(bcs):
        if len(key) != fp_len:
            del bcs[key]
    counts = pd.DataFrame(bcs.most_common())
    return counts.rename(columns={0: 'sequence', 1: 'number'})



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sam_file', help="input sam file")
    parser.add_argument('csv_file', help="output csv file")
    parser.add_argument('--barcode_start')
    parser.add_argument('--barcode_len')
    parser.add_argument('--barcode_id')

    args = parser.parse_args()

    import pandas as pd

    fp_start = int(args.barcode_start)
    fp_len = int(args.barcode_len)

    if not os.path.isfile(args.sam_file):
        print "{} is not a file".format(args.sam_file)
        return 1

    barcodes = None
    if args.barcode_id:
        if not os.path.isfile(args.barcode_id):
            print "{} is not a file".format(args.sam_file)
            return 1
        else:
            barcodes = pd.read_csv(args.barcode_id, header=None, skiprows=[0])
            barcodes = barcodes.rename(columns={0: 'id', 1: 'seq'})

    counts = get_counts(args.sam_file, fp_len, fp_start)

    def id_sequence(seq):
        seq = ''.join(reversed(map(lambda x: {'A':'T','T':'A','C':'G','G':'C','N':'N'}[x],seq)))
        rows = barcodes.ix[barcodes.seq == seq]
        if rows.size > 0:
            return rows.id.iloc[0]
        return 'Unknown'

    if barcodes is not None:
        counts['ids'] = counts.sequence.map(id_sequence)

    counts.to_csv(args.csv_file, index=False)


if __name__ == "__main__":
    import sys

    sys.exit(main())
