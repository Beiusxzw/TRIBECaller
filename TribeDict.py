from .utils import *

class ChromDict(dict):
	def __setitem__(self, key, value):
		if type(key) != int:
			return TypeError("Expected an pysam chromosome identifer, find {}".format(type(key)))
		super().__setitem__(key, value)

class AgDict(dict):
	def __init__(self, chrom, bin_size=100):
		self._chrom = chrom
		self.bin_size = bin_size
	def __setitem__(self, key, value):
		mkey = ROUND_DOWN(key, self.bin_size) if self.bin_size != 1 else key
		if mkey not in self.keys():
			if value == 'G':
				super().__setitem__(mkey, (0, 1, 1))
			elif value == 'A':
				super().__setitem__(mkey, (1, 0, 1))
			else:
				super().__setitem__(mkey, (0, 0, 1))
		else:
			prev_val = self.__getitem__(mkey)
			if value == 'G':
				# mvalue = (prev_val[0]*(prev_val[2]/(prev_val[2]+1)), ((prev_val[1] * prev_val[2] + 1)/(prev_val[2]+1)), prev_val[2]+1)
				mvalue = (prev_val[0], prev_val[1] + 1, prev_val[2]+1)
			elif value == 'A':
				# mvalue = (((prev_val[0] * prev_val[2] + 1)/(prev_val[2]+1)), prev_val[1]*(prev_val[2]/(prev_val[2]+1)), prev_val[2]+1)
				mvalue = (prev_val[0] + 1, prev_val[1], prev_val[2]+1)
			else:
				#mvalue = (prev_val[0]*(prev_val[2]/(prev_val[2]+1)), prev_val[1]*(prev_val[2]/(prev_val[2]+1)), prev_val[2]+1)
				mvalue = (prev_val[0] + 1, prev_val[1], prev_val[2]+1)
			super().__setitem__(mkey, mvalue)

	def merge(self, d2, merge):
	    result = dict()
	    for k,v in self.items():
	        if k in d2:
	            result[k] = merge(v,d2[k])
	    return result