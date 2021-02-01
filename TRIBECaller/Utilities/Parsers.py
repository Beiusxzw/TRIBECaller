# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue

def parse_region(s:str):
	try:
		[chr_,pos] = s.split(":")
		if '-' not in pos:
			return chr_,int(''.join(pos.split(','))),None
		[start,end]=pos.split("-")
		return chr_,int(''.join(start.split(','))),int(''.join(end.split(',')))
	except ValueError:
		print("Invalid genome position. Example: 1:12,000-12,390")

def parse_gtf_line(s):
	chrom, source, feature, start, end, score, strand, frame, attribute=s.split('\t')
<<<<<<< HEAD
	return [chrom, source, feature, start, end, score, strand, frame] + [{x.strip().split(" ")[0]:x.strip().split(" ")[1][1:-1] for x in attribute.split(";")[:-1]}]                 
=======
	return [chrom, source, feature, start, end, score, strand, frame] + [{x.strip().split(" ")[0]:x.strip().split(" ")[1][1:-1] for x in attribute.split(";")[:-1]}]
>>>>>>> commit bug fixed
