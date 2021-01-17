# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Description: Meta classes for TribeCaller
# ------------------------- #
# own Python Modules
# ------------------------- #

from TRIBECaller.Utilities.utils import *
from TRIBECaller.Utilities.PysamExtensions import *

class ChromDict(dict):
	def __setitem__(self, key, value):
		if type(key) != int:
			return TypeError("Expected an pysam chromosome identifer, find {}".format(type(key)))
		super().__setitem__(key, value)

class TribeDict(dict):
	def __init__(self, chrom, bin_size):
		self._chrom = chrom
		self.bin_size = bin_size

	def merge(self, d2, merge):
	    result = dict()
	    for k,v in self.items():
	        if k in d2:
	            result[k] = merge(v,d2[k])
	    return result

	def __call__(self, reads):
		for k,v in zip(reads[0],reads[1]):
			self.__setitem__(k,v)

class NucleotideDict(TribeDict):
	"""
	class MucleotideDict: Extented dictionary for mapping read position and the nucleotides. Debugging use only
	"""
	def __setitem__(self, key, value):
		mkey = ROUND_DOWN(key, self.bin_size) if self.bin_size != 1 else key
		if mkey not in self.keys():
			super().__setitem__(mkey, (value, 1))
		else:
			prev_val = self.__getitem__(mkey)
			super().__setitem__(mkey, (prev_val[0] + value, prev_val[1]+1))

class ASDict(TribeDict):
	"""
	class ASDict: Extented dictionary for mapping A-C/G count and total reads count.
	"""
	def __setitem__(self, key, value):
		mkey = ROUND_DOWN(key, self.bin_size) if self.bin_size != 1 else key
		if mkey not in self.keys():
			if value == 'C' or value == 'G':
				super().__setitem__(mkey, (0, 1, 1))
			elif value == 'A':
				super().__setitem__(mkey, (1, 0, 1))
			else:
				super().__setitem__(mkey, (0, 0, 1))
		else:
			prev_val = self.__getitem__(mkey)
			if value == 'C' or value == 'G':
				# mvalue = (prev_val[0]*(prev_val[2]/(prev_val[2]+1)), ((prev_val[1] * prev_val[2] + 1)/(prev_val[2]+1)), prev_val[2]+1)
				mvalue = (prev_val[0], prev_val[1] + 1, prev_val[2]+1)
			elif value == 'A':
				# mvalue = (((prev_val[0] * prev_val[2] + 1)/(prev_val[2]+1)), prev_val[1]*(prev_val[2]/(prev_val[2]+1)), prev_val[2]+1)
				mvalue = (prev_val[0] + 1, prev_val[1], prev_val[2]+1)
			else:
				#mvalue = (prev_val[0]*(prev_val[2]/(prev_val[2]+1)), prev_val[1]*(prev_val[2]/(prev_val[2]+1)), prev_val[2]+1)
				mvalue = (prev_val[0], prev_val[1], prev_val[2]+1)
			super().__setitem__(mkey, mvalue)

class ATCGDict(TribeDict):
	"""
	class ATCGDict: Extented dictionary for mapping A-T-C-G count and total reads count.
	"""

	def __setitem__(self, key, value):
		mkey = ROUND_DOWN(key, self.bin_size) if self.bin_size != 1 else key
		if mkey not in self.keys():
			if value == 'A':
				super().__setitem__(mkey, (1, 0, 0, 0, 1))
			elif value == 'T':
				super().__setitem__(mkey, (0, 1, 0, 0, 1))
			elif value == 'C':
				super().__setitem__(mkey, (0, 0, 1, 0, 1))
			elif value == 'G':
				super().__setitem__(mkey, (0, 0, 0, 1, 1))
			else:
				super().__setitem__(mkey, (0, 0, 0, 0, 1))
		else:
			prev_val = self.__getitem__(mkey)
			if value == 'A':
				mvalue = (prev_val[0]+1, prev_val[1], prev_val[2], prev_val[3], prev_val[4]+1)
			elif value == 'T':
				mvalue = (prev_val[0], prev_val[1]+1, prev_val[2], prev_val[3], prev_val[4]+1)
			elif value == 'C':
				mvalue = (prev_val[0], prev_val[1], prev_val[2]+1, prev_val[3], prev_val[4]+1)
			elif value == 'G':
				mvalue = (prev_val[0], prev_val[1], prev_val[2], prev_val[3]+1, prev_val[4]+1)
			else:
				mvalue = (prev_val[0], prev_val[1], prev_val[2], prev_val[3], prev_val[4]+1)
			super().__setitem__(mkey, mvalue)
