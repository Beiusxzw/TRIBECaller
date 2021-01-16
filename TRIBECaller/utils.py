# -*- coding: utf-8 -*-

import time
from collections import defaultdict
import pysam

SAM_FLAGS = {"F":[128,80,163,83],"R":[144,64,147,99],"FP":[64,99,80,83],"SP":[144,147,128,163]}

CIGAR = {
0:'M', #BAM_CMATCH
1:'I', #BAM_CINS
2:'D', #BAM_CDEL
3:'N', #BAM_CREF_SKIP
4:'S', #BAM_CSOFT_CLIP
5:'H', #BAM_CHARD_CLIP
6:'P', #BAM_CPAD
7:'=', #BAM_CEQUAL
8:'X', #BAM_CDIFF
9:'B'  #BAM_CBACK
}


def IS_PAIRED(r:pysam.libcalignedsegment.AlignedSegment):
    return r.flag in SAM_FLAGS["F"] or  r.flag in SAM_FLAGS["R"]

def IS_PAIRED_FORWARD(r:pysam.libcalignedsegment.AlignedSegment):
    return r.flag in SAM_FLAGS["F"]

def IS_PAIRED_REVERSE(r:pysam.libcalignedsegment.AlignedSegment):
    return r.flag in SAM_FLAGS["R"]

def IS_PAIRED_FIRST(r:pysam.libcalignedsegment.AlignedSegment):
    return r.flag in SAM_FLAGS["FP"]

def IS_PAIRED_SECOND(r:pysam.libcalignedsegment.AlignedSegment):
    return r.flag in SAM_FLAGS["SP"]

def COMPLEMENT(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def REVERSE_COMPLEMENT(s):
    return COMPLEMENT(s[::-1])

def REVERSE(s):
    return s[::-1]

# IMPLEMENT
def READ_PAIR(bam, reference_id:str, start:int, end:int):
    """
    Generate read without paired information in a BAM file or within a region.
    @args bam: an pysam.AlignmentFile from pysam module
    @args reference_id: should be a string of the contig
    @args start: start position of the fetching
    @args end: end position of the fetching
    @returns a generator of all reads
    """
    pass


def READ_PAIR_GENERATOR(bam, reference_id, start:int, end:int):
    """
    Generate read pairs in a BAM file or within a region,
    Reads are added to read_dict until a pair is found.
    @args bam: an pysam.AlignmentFile from pysam module
    @args reference_id: should be a string of the contig
    @args start: start position of the fetching
    @args end: end position of the fetching
    @returns a generator of tuple containing paired reads
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

def GET_WITHOUT_INSERTION(read):
    forward_sequence = read.get_forward_sequence()
    if 'I' in read.cigarstring:
        s = ''
        ind = 0
        for i in read.cigartuples:
            if i[0] == 0:
                s = s + forward_sequence[ind:ind+i[1]]
                ind += i[1]
            else:
                ind += i[1]
        return s
    else:
        return forward_sequence

def GET_FORWARD_SEQUENCE(read):
    """
    @args an pysam.libcalignedsegment.AlignedSegment object
    """
    if IS_PAIRED_FORWARD(read):
        if IS_PAIRED_FIRST(read):
            return REVERSE_COMPLEMENT(GET_WITHOUT_INSERTION(read))
        else:
            return GET_WITHOUT_INSERTION(read)
    else:
        if IS_PAIRED_FIRST(read):
            return GET_WITHOUT_INSERTION(read)
        else:
            return REVERSE_COMPLEMENT(GET_WITHOUT_INSERTION(read))


def GET_CUR_TIME(s=None):
    return "[{}]".format(time.asctime()) + '\t' + s if s else "[{}]".format(time.asctime()) 

def ROUND_DOWN(a, n):
    return a - a % n

def FLATTEN(x): return [i for s in x for i in s]

def ROUND_UP(a, n):
    return ROUND_DOWN(a + n - 1, n)

def PRINT_LOGO():
    print("""
 _____  __   _____  ___    __  ___      _ _           
/__   \\/__\\  \\_   \\/ __\\  /__\\/ __\\__ _| | | ___ _ __ 
  / /\\/ \\//   / /\\/__\\// /_\\ / /  / _` | | |/ _ \\ '__|
 / / / _  \\/\\/ /_/ \\/  \\//__/ /__| (_| | | |  __/ |   
 \\/  \\/ \\_/\\____/\\_____/\\__/\\____/\\__,_|_|_|\\___|_|                                               
    """)

def PRINT_INFO():
    print("Version: " + "Demo Version")
