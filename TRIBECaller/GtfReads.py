# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#
# ------------------------- #
# own Python Modules
# ------------------------- #

from TRIBECaller.Utilities.Parsers import parse_region, parse_gtf_line

# ------------------------- #
# Python Modules
# ------------------------- #


class GtfReads(object):
	"""
	class GtfReads: Reads a gtf file
	"""
	class Flag:
		"""
		class GtfReads.Flag: Columns name of the gtf file
		"""
		Chromosome = 0
		Source = 1
		Feature = 2
		Start = 3
		End = 4
		Score = 5
		Strand = 6
		Frame = 7
		Attribute =8

	def __init__(self, file_path):
		self.fn = file_path

	def read(self):
		result = []
		with open(self.fn) as f:
			while True: 
				s = f.readline() 
				if not s: 
					break	
				if s[0]=="#": 
					continue 
				s = parse_gtf_line(s)
				yield (s[:8],s[8])

	def fetch(self,region:str=None,gene_ensembl_id:str=None,gene_symbol:str=None,source:str=None):
		"""
		Read gtf files and returns rows contains the provided information.
		@args region: A string describing the genomic region. For example: 1:1280857-1398034
		@args gene_ensembl_id: Ensembl ID of the gene. Optional, if provided, the function will return information of that gene
		@args gene_symbol: Gene Symbol of the gene. Optional,  if provided, the function will return information of that gene
		@args source: returns information with specific source. For example: havana or ensembl
		@returns a list containing the gtf information
		"""
		if region:
			chr_,start,end = parse_region(region)
		result = []
		with open(self.fn) as f:
			while True: 
				s = f.readline() 
				if not s: 
					break	
				if s[0]=="#": 
					continue 
				s = parse_gtf_line(s)
				if region:
					if s[self.Flag.Chromosome] == chr_ and ((int(s[self.Flag.Start]) > start and int(s[self.Flag.Start]) < end) or (int(s[self.Flag.End]) > start and int(s[self.Flag.End]) < end)):
						result.append(s)
				elif gene_ensembl_id:
					if type(gene_ensembl_id) == list:
						if s[self.Flag.Attribute]['gene_id'] in gene_ensembl_id:
							result.append(s)
					else:
						if s[self.Flag.Attribute]['gene_id'] == gene_ensembl_id:
							result.append(s)
				elif gene_symbol:
					if type(gene_symbol) == list:
						if s[self.Flag.Attribute]['gene_name'] in gene_symbol:
							result.append(s)
					else:
						if s[self.Flag.Attribute]['gene_name'] == gene_symbol:
							result.append(s)
				else:
					raise ValueError("")
		return result

