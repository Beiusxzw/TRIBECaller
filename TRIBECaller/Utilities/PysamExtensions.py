# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#
# ------------------------- #
# Python Modules
# ------------------------- #

from collections import defaultdict
from pysam import AlignmentFile, AlignedSegment
from TRIBECaller.Utilities.utils import CLOSEST_INDEX

# Bit-wise flag in sam file format (if set 1):
# Bit Description
# 1   0x1 0b1 template having multiple segments in sequencing
# 2   0x2 0b10    each segment properly aligned according to the aligner
# 4   0x4 0b100   segment unmapped
# 8   0x8 0b1000  next segment in the template unmapped
# 16  0x10    0b10000 SEQ being reverse complemented
# 32  0x20    0b100000    SEQ of the next segment in the template being reverse complemented
# 64  0x40    0b1000000   the first segment in the template
# 128 0x80    0b10000000  the last segment in the template
# 256 0x100   0b100000000 secondary alignment
# 512 0x200   0b1000000000    not passing filters, such as platform/vendor quality controls
# 1024    0x400   0b10000000000   PCR or optical duplicate
# 2048    0x800   0b100000000000  supplementary alignment
# --------------------------------------------------------------------------------------------------------------------------
#      Bits     |  0x1   |     0x2     |    0x4   |       0x8     |   0x10  |     0x20     |      0x40     |       0x80     |
#  Description  | Paired | Proper Pair | unMapped | Mate unMapped | Reverse | Mate Reverse | First in pair | Second in pair | 
# --------------------------------------------------------------------------------------------------------------------------
# Not fequently used
# -------------------------------------------------------------------------------------------------------------------------------
#          0x100        |                  0x200                    |               0x400              |         0x800           |
# Not primary alignment | read fails platform/vendor quality checks | read is PCR or optical duplicate | supplementary alignment |
# -------------------------------------------------------------------------------------------------------------------------------
#
#
# Forward strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
# Reverse strand.
#
# 1. alignments of the first in pair if they map to the forward strand
# 2. alignments of the second in pair if they map to the reverse strand
#
# F: Read Forward
# R: Read Reverse
# FP: Read First in Pair, if paired
# SP: Read Second in Pair, if paired
#
SAM_FLAGS = {"F": [0b10000000,0b01010000,0b10100011,0b01010011,0b00000000],
             "R": [0b10010000,0b01000000,0b10010011,0b01100011,0b00010000],
             "FP":[0b01000000,0b01100011,0b01010000,0b01010011],
             "SP":[0b10010000,0b10010011,0b10000000,0b10100011]}
# Which is equal to 
# {'F': [128, 80, 163, 83, 0],
# 'R' : [144, 64, 147, 99, 16],
# 'FP': [64, 99, 80, 83],
# 'SP': [144, 147, 128, 163]}

"""
Cigar strings indicate the alignment statistics. For full exlanations on CIGAR strings, see: https://drive5.com/usearch/manual/cigar.html
"""
CIGAR = {
0:'M', #BAM_CMATCH; Match (alignment column containing two letters). This could contain two different letters (mismatch) or two identical letters. USEARCH generates CIGAR strings containing Ms rather than X's and ='s (see below).
1:'I', #BAM_CINS; Insertion (gap in the query sequence).
2:'D', #BAM_CDEL; Deletion (gap in the target sequence).
3:'N', #BAM_CREF_SKIP; 
4:'S', #BAM_CSOFT_CLIP; Segment of the query sequence that does not appear in the alignment. This is used with soft clipping, where the full-length query sequence is given (field 10 in the SAM record). In this case, S operations specify segments at the start and/or end of the query that do not appear in a local alignment.
5:'H', #BAM_CHARD_CLIP; Segment of the query sequence that does not appear in the alignment. This is used with hard clipping, where only the aligned segment of the query sequences is given (field 10 in the SAM record). In this case, H operations specify segments at the start and/or end of the query that do not appear in the SAM record.
6:'P', #BAM_CPAD
7:'=', #BAM_CEQUAL; Alignment column containing two identical letters. USEARCH can read CIGAR strings using this operation, but does not generate them.
8:'X', #BAM_CDIFF; Alignment column containing a mismatch, i.e. two different letters. USEARCH can read CIGAR strings using this operation, but does not generate them.
9:'B'  #BAM_CBACK
}


def MERGE_DICT(d1,d2,merge):
    """
    Merge two dictionaries by key
    """
    result = dict()
    for k,v in d1.items():
        if k in d2:
            result[k] = merge(v,d2[k])
        else:
            result[k] = v
    return result

def MERGE_AS_DICT(d1,d2):
    """
    Merge two dictionaries by key
    """
    result = dict()
    for k,v in d1.items():
        if k in d2:
            result[k] = (v[0]+d2[k][0],v[1]+d2[k][1],v[2]+d2[k][2])
        else:
            result[k] = v
    return result

def IS_PAIRED(r:AlignedSegment):
    """
    return True if the AlignedSegment is a paired-read
    """
    return r.flag & 0b10000000 or r.flag & 0b01000000

def IS_EQ_REF_FORWARD(r:AlignedSegment):
    """
    return True if the read of AlignedSegment is the same as the forward strand
    Don't confused with read-forward in the bit-wise flag
    """
    return r.flag in SAM_FLAGS["F"]

def IS_EQ_REF_REVERSE(r:AlignedSegment):
    """
    return True if the read of AlignedSegment is is the same as the reverse strand
    Don't confused with read-forward in the bit-wise flag
    """
    return r.flag in SAM_FLAGS["R"]

def IS_PAIRED_FIRST(r:AlignedSegment):
    """
    return True if the AlignedSegment is Read-First-in-Pair
    return False if the AlignedSegment is not Read-First-in-Pair or the read do not come from pair-ended library
    """
    return r.flag & 0b01000000
    # return r.flag in SAM_FLAGS["FP"]

def IS_PAIRED_SECOND(r:AlignedSegment):
    """
    return True if the AlignedSegment is Read-Second-in-Pair
    return False if the AlignedSegment is not Read-Second-in-Pair or the read do not come from pair-ended library
    """
    return r.flag & 0b10000000
    # return r.flag in SAM_FLAGS["SP"]

def COMPLEMENT(seq):
    """
    return Complement sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def REVERSE_COMPLEMENT(s):
    """
    return Reverse Complement sequence
    """
    return COMPLEMENT(s[::-1])

def REVERSE(s):
    return s[::-1]

def READ_PAIR(read1, read2, keep=1):
    assert(len(read1) == len(read2))
    if read1[0][0] <= read2[0][0]:
        if read1[0][-1] >= read2[0][0]:
            # ----------------------> R1
            #              <---------------------R 2
            if keep == 1:
                try:
                    index = read2[0].index(read1[0][-1])
                except:
                    index = CLOSEST_INDEX(read2[0], read1[0][-1])
                return (read1, (read2[0][ index :], read2[1][ index :]))
            else:
                try:
                    index = read1[0].index(read2[0][0])
                except:
                    index = CLOSEST_INDEX(read1[0], read2[0][0])
                return ((read1[0][ : index ], read1[1][: index ]), read2)

        else:
            # ----------------------> R1
            #                        <---------------------R 2
            return (read1,read2)
    else:
        if read2[0][-1] >= read1[0][0]:
            #               ----------------------> R1
            # <---------------------R 2
            if keep == 1:
                try:
                    index = read2[0].index(read1[0][0])
                except:
                    index = CLOSEST_INDEX(read2[0], read1[0][0])

                return (read1, (read2[0][: index], read2[1][:index]))
            else:
                try:
                    index = read1[0].index(read1[0][-1])
                except:
                    index = CLOSEST_INDEX(read1[0], read1[0][-1])
                return ((read1[0][index:], read1[1][index:]), read2)
        else:
            #                        ----------------------> R1
            # <---------------------R 2            
            return (read1,read2)

def READ_PAIR_GENERATOR(bam, reference_id, start:int, end:int):
    """
    Generate read pairs in a BAM file or within a region,
    Reads are added to read_dict until a pair is found.
    @args bam: an pysam.AlignmentFile from pysam module
    @args reference_id: should be a string of the contig
    @args start: start position of the fetching
    @args end: end position of the fetching
    @returns a generator object of tuple containing paired reads
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(reference_id,start,end):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def GET_WITHOUT_INSERTION(read, rev, clip=False,debug=False):
    forward_sequence = read.get_forward_sequence()
    s = ''
    ind = 0
    ct = read.cigartuples if rev else read.cigartuples[::-1]
    for i in ct:
        if i[0] == 0:
            s = s + forward_sequence[ind:ind+i[1]]
            ind += i[1]
        elif i[0] == 1:
            ind += i[1]
        elif i[0] == 2:
            s = s + forward_sequence[ind:ind+i[1]]
            ind += i[1]
        elif i[0] == 4 and clip:
            ind += i[1]
        else:
            if debug:
                raise ValueError("mismatch other than Deletion or Insertion are not accepted! Please check your alignment")
            else:
                return None
    return s

def GET_FORWARD_SEQUENCE(read):
    """
    @args an AlignedSegment object
    """
    if not GET_WITHOUT_INSERTION(read, rev = False):
        return None
    if IS_PAIRED(read):
        if IS_EQ_REF_FORWARD(read):
            if IS_PAIRED_FIRST(read):
                return REVERSE_COMPLEMENT(GET_WITHOUT_INSERTION(read, rev = False))
            else:
                return GET_WITHOUT_INSERTION(read, rev=True)
        else:
            if IS_PAIRED_FIRST(read):
                return GET_WITHOUT_INSERTION(read,rev=True)
            else:
                return REVERSE_COMPLEMENT(GET_WITHOUT_INSERTION(read,rev=False))
    else:
        if IS_EQ_REF_FORWARD(read):
            return GET_WITHOUT_INSERTION(read, rev=False)
        else:
            return REVERSE_COMPLEMENT(GET_WITHOUT_INSERTION(read, rev=True))