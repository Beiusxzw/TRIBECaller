from .TribeDict import *
import pysam

class TribeReads(object):
	"""
	class description: A wrapper for AlignmentFile from pysam module that provides computation
	                   of the A-G content
	@author: Ziwei Xue
	"""
	def __init__(self, file_path:str=None, pysam_file = None):
		"""
		@args file_path: A string of the input file path. Should be suffixed with bam. If no index 
		                 has been built, pysam will print a warning message.
		@args pysam_file: An AlignmentFile Object can be passed into the wrapper if no file path is
		                  provided
		"""
		super(TribeReads, self).__init__()
		self.aligned_reads = pysam.AlignmentFile(file_path) if not pysam_file else pysam_file
		self._header = self.aligned_reads.header.as_dict()
		self._chrom_sizes = dict(list(map(lambda x:list(x.values()), self._header["SQ"])))

	def build_ag_dict(self, reference_id, start:int, end:int, bin_size = 100):
		"""
		@args reference_id: The reference id of the contig. Can be a integer flag from the pysam or
		                    the real chromosome string, such as "chr1"
		@args start: the start of the position of the chromosome where reads will be mapped
		@args end: the end of the position of the chromosome where reads will be mapped
		@args bin_size: window of the computed G-C content 
		@returns An AgDict that contains information of the a-g content
		"""
		if type(reference_id) == int:
			ref_name = self.aligned_reads.get_reference_name(reference_id)
		else:
			ref_name = reference_id
		ag_dict = AgDict(ref_name, bin_size)
		for read in self.aligned_reads.fetch(ref_name, start, end):
			for pos,nuc in zip(read.get_reference_positions(), read.get_forward_sequence()):
				ag_dict[pos] = nuc
		return ag_dict

	async def build_ag_dict_async(self, reference_id, start:int, end:int, bin_size = 100):
		"""
		An asynchronous function of build_ag_dict
		"""
		if type(reference_id) == int:
			ref_name = self.aligned_reads.get_reference_name(reference_id)
		else:
			ref_name = reference_id
		ag_dict = AgDict(ref_name, bin_size)
		for read in self.aligned_reads.fetch(ref_name, start, end):
			for pos,nuc in zip(read.get_reference_positions(), read.get_forward_sequence()):
				ag_dict[pos] = nuc
		return ag_dict

	
	def get_chrom_size(self, reference_id):
		if type(reference_id) == int:
			return self._chrom_sizes[self.aligned_reads.get_reference_name(reference_id)]
		else:
			return self._chrom_sizes[reference_id]

	def get_chrom_name():
		return self._chrom_sizes.kyes()

