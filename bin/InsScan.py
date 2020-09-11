import pysam
import sys
import argparse

class Options:
    def __init__(self):
        self.parser = argparse.ArgumentParser("Detect short and long insertion")
        self.parser.add_argument('-bam',help='bam file',dest='bam')
        self.parser.add_argument('-n',help='sample name',dest='name')
        self.parser.add_argument('-r',help='region <chr:start-end>',dest='region')
        self.parser.add_argument('-outdir',help='output dir',dest='outdir')


        args = self.parser.parse_args()

        if args.bam:
            self.bam = args.bam
        else:
            self.parser.error('BAM file not specified')

        if args.name:
            self.name = args.name
        else:
            self.parser.error('sample name not specified')

        if args.outdir:
            self.outdir = args.outdir
        else:
            self.parser.error('outdir not specified')


def main():
    options = Options()
    samfile = pysam.AlignmentFile(options.bam,"rb")
    for read in samfile.fetch():
        print(read)
        # return an aligned read
        #if read.is_duplicate or read.qcfail:
            #next
        # bwa mem's -M will trans supplemantary align into secondary. sup/secondary read's CIGAR will contain soft-clip or hard-clip info, these reads will be used to infer SV.

        # if this read's CIGAR contain 'I'
        #if read.cigarstring





if __name == "__main__":
    main()



