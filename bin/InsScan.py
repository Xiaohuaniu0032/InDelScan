import pysam
import sys
import argparse
import re

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


cigar_value = ['M','I','D','N','S','H','P','=','X']

def main():
    options = Options()
    samfile = pysam.AlignmentFile(options.bam,"rb")

    outfile = options.outdir + '/%s.softclip.read.txt' % (options.name)
    of = open(outfile,'w')
    header = "seqname\tchr\tleft_align_pos\tcigar\tinfer_left_breakpoint\tinfer_right_breakpoint\tseq\tmapq\tflags\tstrand\tread1/2\n";
    of.write(header)

    for read in samfile.fetch():
        if read.is_duplicate or read.is_qcfail or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < 20:
            continue

        print(read)
        print(read.reference_name)
        print(read.reference_start)
        print(read.cigarstring)
        print(read.cigartuples)
        #print(read.query_alignment_end)
        #print(read.query_alignment_length)
        #print(read.get_cigar_stats())
        #print(read.infer_query_length())
        #print(read.infer_read_length())

        cigar = read.cigarstring
        chrom = read.reference_name
        left_aln_pos = read.reference_start + 1
        seq = read.query_sequence
        mapq = read.mapping_quality
        baseQ = read.query_qualities
        seqname = read.query_name
        flags = read.flag
        cigar_string = read.cigarstring
        cigar_tuple = read.cigartuples
        
        '''
        * return an aligned read
        * bwa mem's -M will trans supplemantary align into secondary. sup/secondary read's CIGAR will contain soft-clip or hard-clip info, these reads will be used to infer SV.
        '''

        # strand
        if read.is_reverse:
            strand = '-'
        else:
            strand = '+'

        # r1 or r2
        if read.is_read1:
            r1r2 = "read1"
        else:
            r1r2 = "read2"

        # skip read with (\d+)S(\d+)M(\d+)S
        # count S num
        s_num = 0
        for i in cigar_tuple:
            if i[0] == 4:
                s_num += 1

        if s_num == 2:
            # can print some log info
            continue


        # check if exists I in cigar
        ins_num = 0
        for i in cigar_tuple:
            if i[0] == 1:
                ins_num += 1

        #print("soft-clip and ins info is:")
        #print(s_num)
        #print(ins_num)
        # skip read without S or I for that these reads can not used to infer short/long Insert
        if s_num == 0 and ins_num == 0:
            continue

        break_point = infer_breakpoint(read)

        val = seqname + '\t' + str(chrom) + '\t' + str(left_aln_pos) + '\t' + cigar_string + '\t' + str(break_point[0]) + '\t' + str(break_point[1]) + '\t' + seq + '\t' + str(mapq) + '\t' + str(flags) + '\t' + strand + '\t' + r1r2

        of.write(val+'\n')
    of.close()

def infer_breakpoint(readAlignedSegment):
    '''
    M   BAM_CMATCH      0
    I   BAM_CINS        1
    D   BAM_CDEL        2
    N   BAM_CREF_SKIP   3
    S   BAM_CSOFT_CLIP  4
    H   BAM_CHARD_CLIP  5
    P   BAM_CPAD        6
    =   BAM_CEQUAL      7
    X   BAM_CDIFF       8
    B   BAM_CBACK       9
    
    M/D/N/=/X will consume ref

    returned value is the consumed ref length
    '''

    #CIGAR_FORMAT = re.compile(r'(\d+)M(\d+)I(\d+)M')
    #    if CIGAR_FORMAT.match(read.cigarstring):


    left_aln_pos = readAlignedSegment.reference_start + 1
    cigar_tuple = readAlignedSegment.cigartuples
    cigar_string = readAlignedSegment.cigarstring
    #print(cigar_string)

    
    # check if exists I in cigar
    ins_num = 0
    for i in cigar_tuple:
        if i[0] == 1:
            ins_num += 1

    if ins_num == 0:
        # no ins
        # left or right
        s_idx = cigar_string.index('S')
        if s_idx + 1 == len(cigar_string):
            # right
            ref_consume_len = 0
            for i in cigar_tuple:
                if i[0] == 0 or i[0] == 2 or i[0] == 3 or i[0] == 7 or i[0] == 8:
                    ref_consume_len += i[1]
            left_break_point = left_aln_pos + ref_consume_len - 1
            right_break_point = 'NA'
        else:
            # left
            left_break_point = 'NA'
            right_break_point = left_aln_pos
    else:
        # has ins
        # which pos I is
        flag = 0
        for i in cigar_tuple:
            flag += 1
            if i[0] == 1:
                I_pos = flag
                #print('I pos is:')
                #print(I_pos)
                break

        I_pos = flag
        flag = 0
        ref_consume_len = 0
        for i in cigar_tuple:
            flag += 1
            if flag <= I_pos:
                # how many ref bp consumed before Ins pos
                if i[0] == 0 or i[0] == 2 or i[0] == 3 or i[0] == 7 or i[0] == 8:
                    ref_consume_len += i[1]
                    #print("ref consumed len is:")
                    #print(ref_consume_len)

        left_break_point = left_aln_pos + ref_consume_len - 1
        right_break_point = left_break_point + 1


    val = (left_break_point,right_break_point)
    return(val)
    

if __name__ == "__main__":
    main()



