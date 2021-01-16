# -*- coding: utf-8 -*-

from TRIBECaller.TribeEditingCaller import *
from TRIBECaller.utils import *
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-t', '--target', type=str, help='input taregt bam file')
parser.add_argument('-c', '--control', type=str, help='input control bam file')
parser.add_argument('-f', '--frameSize', type=int, default=600000,help='bin size of calling editing event')
parser.add_argument('-b', '--binSize', type=int, default=1,help='bin size of calling editing event')
parser.add_argument('-o', '--outPrefix', type=str, default = "EditingEvents", help='prefix of output bed file')
args = parser.parse_args()

def call_editing_sites():
    PRINT_LOGO()
    PRINT_INFO()
    TEC = TribeCaller(args.target,args.control)
    result = TEC()
    res = []
    with open(args.outPrefix + "bed", "w+") as f:
        for k,v in result.items():
            for j in v:
                f.write("\t".join([k,j[0]+1,j[0]+2,j[1],j[2]]))