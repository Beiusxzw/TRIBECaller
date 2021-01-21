# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#
# ------------------------- #
# own Python Modules
# ------------------------- #

from TRIBECaller.TribeDict import *
from TRIBECaller.Utilities.MapReduce import *
# ------------------------- #
# Python Modules
# ------------------------- #

import pysam
from functools import partial
import time

class TribeReads(object):
	"""
	class TribeReads: A wrapper for pysam.AlignmentFile from pysam module that provides computation
					   of the A-G content
	"""
	def __init__(self, file_path:str=None, executor=None, chunk_size=10000, num_workers=12):
		"""
		@args file_path: A string of the input file path. Should be suffixed with bam. If no index 
						 has been built, pysam will print a warning message.
		"""
		super(TribeReads, self).__init__()
		self.file_path = file_path
		self._header = pysam.AlignmentFile(self.file_path).header.as_dict()
		self._chrom_sizes = dict(list(map(lambda x:list(x.values()), self._header["SQ"])))
		self.executor = executor
		self.num_workers = num_workers
		self.chunk_size = chunk_size

	def get_reference_id(self, reference_id):
		if type(reference_id) == int:
			return pysam.AlignmentFile(self.file_path).get_reference_name(reference_id)
		else:
			return reference_id

	def fetch(self,*args,**kwargs):
		return pysam.AlignmentFile(self.file_path).fetch(*args,**kwargs)

	def get_reads(self, reference_id, start:int, end:int):
		reference_id = self.get_reference_id(reference_id)
		return list(filter(lambda x:x[1], map(lambda read: (read.get_reference_positions(), GET_FORWARD_SEQUENCE(read)), pysam.AlignmentFile(self.file_path).fetch(reference_id, start, end))))

	def get_reads_paired(self, reference_id, start:int, end:int, exclude=False):
		"""
		@args reference_id
		@args start
		@args end
		@args exclude
		"""
		reference_id = self.get_reference_id(reference_id)
		return FLATTEN(list(map(lambda reads: [(read.get_reference_positions(), GET_FORWARD_SEQUENCE(read)) for read in reads], READ_PAIR_GENERATOR(pysam.AlignmentFile(self.file_path), reference_id, start, end))))

	def map_func(self, chrom, bin_size, reads):
		s = time.time()
		d = ASDictPar(chrom, bin_size)
		for i in reads:
			d(i)
		print("Mapped :{}".format(time.time()-s))
		return d

	def reduce_func(self, ds):
		s = time.time()
		if not ds:
			return
		if len(ds) == 1:
			print("Reduced :{}".format(time.time()-s))
			return ds[0]
		keys = {key for d in ds for key in d}
		td = TupleDict()
		for key in keys:
			for d in ds:
				if key in d:
					td.__setitem__(key,d[key]) 
		print("Reduced ha :{}".format(time.time()-s))
		return td

	def build_nucleotides_dict(self, reference_id, start:int, end:int, bin_size = 1):
		reference_id = self.get_reference_id(reference_id)
		reads = self.get_reads(reference_id, start, end)
		nuc_dict = NucleotideDict(reference_id, bin_size)
		for i in reads:
			nuc_dict(i)
		return nuc_dict

	def build_nucleotides_dict_par(self, reference_id, start:int, end:int, bin_size = 1):
		raise NotImplementedError("This method has not been implemented")

	def build_as_dict(self, reference_id, start:int, end:int, bin_size = 1):
		"""
		Build a ASDict according to reads in a genomic region with certain bin size
		@args reference_id: The reference id of the contig. Can be a integer flag from the pysam or
							the real chromosome string, such as "chr1"
		@args start: the start of the position of the chromosome where reads will be mapped
		@args end: the end of the position of the chromosome where reads will be mapped
		@args bin_size: window of the computed A-(G/C) content 
		@returns An ASDict that contains information of the A-(G/C) content
		"""
		reference_id = self.get_reference_id(reference_id)
		reads = self.get_reads(reference_id, start, end)
		as_dict = ASDict(reference_id, bin_size)
		for i in reads:
			as_dict(i)
		return as_dict

	def build_as_dict_par(self, reference_id, start:int, end:int, bin_size = 1):
		reference_id = self.get_reference_id(reference_id)
		reads = self.get_reads(reference_id, start, end)
		mr = MapReduce(partial(self.map_func,reference_id,bin_size), self.reduce_func, self.chunk_size, self.num_workers)
		return mr.proc_map_reduce(reads)
		
	def build_atcg_dict(self, reference_id, start:int, end:int, bin_size = 1):
		reference_id = self.get_reference_id(reference_id)
		reads = self.get_reads(reference_id, start, end)
		atcg_dict = ATCGDict(reference_id, bin_size)
		for i in reads:
			atcg_dict(i)
		return atcg_dict
	
	def build_atcg_dict_par(self, reference_id, start:int, end:int, bin_size = 1):
		raise NotImplementedError("This method has not been implemented")

	def get_chrom_size(self, reference_id):
		if type(reference_id) == int:
			return self._chrom_sizes[pysam.AlignmentFile(self.file_path).get_reference_name(reference_id)]
		else:
			return self._chrom_sizes[reference_id]

	def get_chrom_name(self):
		return self._chrom_sizes.keys()

