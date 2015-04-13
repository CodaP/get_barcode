import pysam
import pandas as pd
import os
from collections import Counter
 
import argparse
 
 
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sam_file', help="input sam file")
    parser.add_argument('csv_file', help="output csv file")
    parser.add_argument('--barcode_start')
    parser.add_argument('--barcode_len')
 
    args = parser.parse_args()
 
    fp_start = int(args.barcode_start)
    fp_len = int(args.barcode_len)
 
    if not os.path.isfile(args.sam_file):
        print "{} is not a file".format(args.sam_file)
        return 1
 
    sam = pysam.Samfile(args.sam_file)
 
    bcs = Counter([x.seq[(fp_start - x.reference_start - 1):][:fp_len] for x in sam if x.reference_start <= fp_start and x.reference_end >= fp_start+fp_len])
 
    for key in list(bcs):
        if len(key) != fp_len:
            del bcs[key]
 
    counts = pd.DataFrame(bcs.most_common())
 
    counts.to_csv(args.csv_file, index=False, header=['sequence', 'number'])
 
if __name__ == "__main__":
    import sys
    sys.exit(main())
