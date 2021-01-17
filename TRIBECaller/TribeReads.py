# -*- coding: utf-8 -*-
#!/usr/bin/env python
# ------------------------- #
# own Python Modules
# ------------------------- #

from TRIBECaller.TribeDict import *

# ------------------------- #
# Python Modules
# ------------------------- #

import pysam
from functools import partial
from concurrent.futures import ProcessPoolExecutor   

class TribeReads(object):
	"""
	class TribeReads: A wrapper for pysam.AlignmentFile from pysam module that provides computation
					   of the A-G content
	"""
	def __init__(self, file_path:str=None,  executor=None):
		"""
		@args file_path: A string of the input file path. Should be suffixed with bam. If no index 
						 has been built, pysam will print a warning message.
		"""
		super(TribeReads, self).__init__()
		self.file_path = file_path
		self._header = pysam.AlignmentFile(self.file_path).header.as_dict()
		self._chrom_sizes = dict(list(map(lambda x:list(x.values()), self._header["SQ"])))
		self.executor = executor

	def fetch(self,*args,**kwargs):
		return pysam.AlignmentFile(self.file_path).fetch(*args,**kwargs)

	def get_reads(self, reference_id, start:int, end:int):
		if type(reference_id) == int:
			ref_name = pysam.AlignmentFile(self.file_path).get_reference_name(reference_id)
		else:
			ref_name = reference_id
		return list(map(lambda read: (read.get_reference_positions(), GET_FORWARD_SEQUENCE(read)), pysam.AlignmentFile(self.file_path).fetch(ref_name, start, end)))

	def get_reads_paired(self, reference_id, start:int, end:int):
		if type(reference_id) == int:
			ref_name = pysam.AlignmentFile(self.file_path).get_reference_name(reference_id)
		else:
			ref_name = reference_id

		return FLATTEN(list(map(lambda reads: [(read.get_reference_positions(), GET_FORWARD_SEQUENCE(read)) for read in reads], READ_PAIR_GENERATOR(pysam.AlignmentFile(self.file_path), ref_name, start, end))))

	def build_nucleotides_dict(self, reference_id, start:int, end:int, bin_size = 1):
		if type(reference_id) == int:
			ref_name = pysam.AlignmentFile(self.file_path).get_reference_name(reference_id)
		else:
			ref_name = reference_id
		reads = self.get_reads(ref_name, start, end)
		nuc_dict = NucleotideDict(ref_name, bin_size)
		for i in reads:
			nuc_dict(i)
		return nuc_dict

	def build_as_dict(self, reference_id, start:int, end:int, bin_size = 1):
		"""
		@args reference_id: The reference id of the contig. Can be a integer flag from the pysam or
							the real chromosome string, such as "chr1"
		@args start: the start of the position of the chromosome where reads will be mapped
		@args end: the end of the position of the chromosome where reads will be mapped
		@args bin_size: window of the computed G-C content 
		@returns An ASDict that contains information of the t-c content
		"""
		if type(reference_id) == int:
			ref_name = pysam.AlignmentFile(self.file_path).get_reference_name(reference_id)
		else:
			ref_name = reference_id
		reads = self.get_reads(ref_name, start, end)
		as_dict = ASDict(ref_name, bin_size)
		for i in reads:
			as_dict(i)
		return as_dict

	def build_atcg_dict(self, reference_id, start:int, end:int, bin_size = 1):
		if type(reference_id) == int:
			ref_name = pysam.AlignmentFile(self.file_path).get_reference_name(reference_id)
		else:
			ref_name = reference_id
		reads = self.get_reads(ref_name, start, end)
		atcg_dict = ATCGDict(ref_name, bin_size)
		for i in reads:
			atcg_dict(i)
		return atcg_dict
	
	def get_chrom_size(self, reference_id):
		if type(reference_id) == int:
			return self._chrom_sizes[pysam.AlignmentFile(self.file_path).get_reference_name(reference_id)]
		else:
			return self._chrom_sizes[reference_id]

	def get_chrom_name(self):
		return self._chrom_sizes.keys()

