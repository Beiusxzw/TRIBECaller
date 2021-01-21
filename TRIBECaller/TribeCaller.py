# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#
# ------------------------- #
# own Python Modules
# ------------------------- #

from TRIBECaller.TribeReads import *
from TRIBECaller.Utilities.utils import *
from TRIBECaller.Utilities.PysamExtensions import *

# ------------------------- #
# Python Modules
# ------------------------- #

import numpy as np 
import scipy.stats as stats
from typing import Optional, Union, Tuple
import tqdm

# ------------------------- #
# Main Class
# ------------------------- #

class TribeCaller(object):
	def __init__(self, target_reads_path:str, 
				 control_reads_path:str, 
				 frame_size:int = 600000, 
				 bin_size:int = 1,
				 num_threads:int=1,
				 run_async: bool = False):
		super(TribeCaller, self).__init__()
		self.target = TribeReads(target_reads_path)
		self.control = TribeReads(control_reads_path)
		if len(self.target._chrom_sizes) != len(self.control._chrom_sizes):
			print(GET_CUR_TIME("The length of chromosome number is not the same between target and control"))
		for i in self.target._chrom_sizes.keys():
			if self.target._chrom_sizes[i] != self.control._chrom_sizes[i]:
				print(GET_CUR_TIME("Different chromosome size found. Please check if use the same aligner"))
		self._chrom_sizes = self.target._chrom_sizes
		self.frame_size = frame_size
		self.bin_size = bin_size
		self.num_threads = num_threads

	def compute_atcg_content(self, reference_id, start:int, end:int):
		target_atcgdict =  self.target.build_atcg_dict(reference_id, start, end, self.bin_size)
		control_atcgdict = self.control.build_atcg_dict(reference_id, start, end, self.bin_size)
		return target_atcgdict.merge(control_atcgdict, lambda x,y:(x,y))

	def compute_region_as_percentage(self,reference_id,start:int,end:int=None, threshold=2, content_threshold=0.8,pvalue_cutoff=0.05,diff_cutoff=None,call_editing_sites=False):
		if end:
			if end - start > 200:
				raise TypeError("the region length should not exceed 200bp")
		else:
			start = start - 20
			end = start + 40
		target_atcgdict =  self.target.build_atcg_dict(reference_id, start, end, self.bin_size)
		target_atcgdict = {k: target_atcgdict[k] if k in target_atcgdict.keys() else (0,0,0,0,0) for k in range(start,end)}
		control_atcgdict = self.control.build_atcg_dict(reference_id, start, end, self.bin_size)
		control_atcgdict = {k: control_atcgdict[k] if k in control_atcgdict.keys() else (0,0,0,0,0) for  k in range(start,end)}
		count,result = [],[]
		for i in range(start,end):
			t,c = target_atcgdict[i],control_atcgdict[i]
			if call_editing_sites:
				result.append([[t[0]/t[4], t[1]/t[4],t[2]/t[4],t[3]/t[4],t[4]] if t[4] > 0 else [0,0,0,0,0] ,[c[0]/c[4], c[1]/c[4],c[2]/c[4],c[3]/c[4],c[4]] if c[4] > 0 else [0,0,0,0,0]])
				count.append([[t[0], t[1],t[2],t[3],t[4]] if t[4] > 0 else [0,0,0,0,0] ,[c[0], c[1],c[2],c[3],c[4]] if c[4] > 0 else [0,0,0,0,0]])
			else:
				result.append([[t[0]/t[4], t[1]/t[4],t[2]/t[4],t[3]/t[4]] if t[4] > 0 else [0,0,0,0] ,[c[0]/c[4], c[1]/c[4],c[2]/c[4],c[3]/c[4]] if c[4] > 0 else [0,0,0,0]])
		if call_editing_sites:
			res = list(map(lambda x: (stats.fisher_exact([[x[0][3], x[1][3]],[x[0][0], x[1][0]]], "greater"), stats.fisher_exact([[x[0][2], x[1][2]],[x[0][1], x[1][1]]], "greater")), count))
			odds_f, pval_f, odds_r, pval_r = np.array(list(map(lambda x:x[0], map(lambda x:x[0],res)))),np.array(list(map(lambda x:x[1], map(lambda x:x[0],res)))),np.array(list(map(lambda x:x[0], map(lambda x:x[1],res)))),np.array(list(map(lambda x:x[1], map(lambda x:x[1],res))))
			with np.errstate(invalid='ignore'):
				return list(range(start,end)), list(map(lambda x:[x[0][:4],x[1][:4]], result)), ((odds_f > threshold) & (pval_f < pvalue_cutoff) & (np.array(list(map(lambda x:x[1][0]/x[1][4]  > content_threshold if x[1][4] > 0 else False, count))))) | ((odds_r > threshold) & (pval_r < pvalue_cutoff) & (np.array(list(map(lambda x:x[1][1]/x[1][4] > content_threshold if x[1][4] > 0 else False, count)))))
		return list(range(start,end)), result

	def call_editing_region(self,reference_id, start:int, end:int, threshold=2,content_threshold=0.8,pvalue_cutoff=0.05, complex=False):
		merged_atcgdict = self.compute_atcg_content(reference_id, start, end)
		values=list(merged_atcgdict.values())
		diff = list(map(lambda x: ((x[0][3]/(x[0][0]+x[0][3])) - (x[1][3]/(x[1][0]+x[1][3])) if (x[0][0]+x[0][3] > 0 and x[1][0]+x[1][3] > 0) else None, (x[0][2]/(x[0][1]+x[0][2]) - (x[1][2]/(x[1][1]+x[1][2]))) if (x[0][1]+x[0][2] > 0 and x[1][1]+x[1][2] > 0) else None ), values))
		res = list(map(lambda x: (stats.fisher_exact([[x[0][3], x[1][3]],[x[0][0], x[1][0]]], "greater"), stats.fisher_exact([[x[0][2], x[1][2]],[x[0][1], x[1][1]]], "greater")), values))
		odds_f, pval_f, odds_r, pval_r = np.array(list(map(lambda x:x[0], map(lambda x:x[0],res)))),np.array(list(map(lambda x:x[1], map(lambda x:x[0],res)))),np.array(list(map(lambda x:x[0], map(lambda x:x[1],res)))),np.array(list(map(lambda x:x[1], map(lambda x:x[1],res))))
		if complex:
			out = []
			for x in range(len(res)):
				if (odds_f[x] > threshold and pval_f[x] < pvalue_cutoff and values[x][1][0]/values[x][1][4] > content_threshold):
					out.append((list(merged_atcgdict.keys())[x],diff[x][0],values[x][0][0],values[x][0][3],values[x][1][0],values[x][1][3], pval_f[x], '+'))
				if (odds_r[x] > threshold and pval_r[x] < pvalue_cutoff and values[x][1][1]/values[x][1][4] > content_threshold):
					out.append((list(merged_atcgdict.keys())[x],diff[x][1],values[x][0][1],values[x][0][2],values[x][1][1],values[x][1][2], pval_r[x], '-'))
			return out
		else:
			out = []
			for x in range(len(res)):
				if (odds_f[x] > threshold and pval_f[x] < pvalue_cutoff and values[x][1][0]/values[x][1][4] > content_threshold):
					out.append((list(merged_atcgdict.keys())[x],diff[x][0],pval_f[x],'+'))
				if (odds_r[x] > threshold and pval_r[x] < pvalue_cutoff and values[x][1][1]/values[x][1][4] > content_threshold):
					out.append((list(merged_atcgdict.keys())[x],diff[x][1],pval_r[x],'-'))
			return out

	def call_editing_region_coverage(self,reference_id, start:int, end:int, threshold=2, content_threshold=0.8,pvalue_cutoff=0.05):
		def _helper(x):
			return 0 if x < 0 else x
		print(GET_CUR_TIME("Computing bam coverage and nucleotides content"))
		#if end - start > 200000:
		#	raise TypeError("the region length should not exceed 200000bp")
		merged_atcgdict = self.compute_atcg_content(reference_id, start, end)
		values=list(merged_atcgdict.values())
		diff = list(map(lambda x: ( _helper((x[0][3]/(x[0][0]+x[0][3])) - (x[1][3]/(x[1][0]+x[1][3]))) if (x[0][0]+x[0][3] > 0 and x[1][0]+x[1][3] > 0) else 0, _helper((x[0][2]/(x[0][1]+x[0][2]) - (x[1][2]/(x[1][1]+x[1][2])))) if (x[0][1]+x[0][2] > 0 and x[1][1]+x[1][2] > 0) else 0) if x[0][4] > 5 and x[1][4] > 5 else 0, values))
		res = list(map(lambda x: (stats.fisher_exact([[x[0][3], x[1][3]],[x[0][0], x[1][0]]], "greater"), stats.fisher_exact([[x[0][2], x[1][2]],[x[0][1], x[1][1]]], "greater")), values))
		odds_f, pval_f, odds_r, pval_r = np.array(list(map(lambda x:x[0], map(lambda x:x[0],res)))),np.array(list(map(lambda x:x[1], map(lambda x:x[0],res)))),np.array(list(map(lambda x:x[0], map(lambda x:x[1],res)))),np.array(list(map(lambda x:x[1], map(lambda x:x[1],res))))
		with np.errstate(invalid='ignore'):
			return list(merged_atcgdict.keys()), list(map(lambda y: (y[0][4],y[1][4]), merged_atcgdict.values())), ((odds_f > threshold) & (pval_f < pvalue_cutoff) & (np.array(list(map(lambda x:x[1][0]/x[1][4] > content_threshold, values))))) | ((odds_r > threshold) & (pval_r < pvalue_cutoff) & (np.array(list(map(lambda x:x[1][1]/x[1][4] > content_threshold, values))))), list(map(lambda x:max(x) if type(x) == tuple else x, diff))

	def compute_nucleotides_coverage(self,reference_id, start:int, end:int, threshold=2):
		merged_atcgdict = self.compute_atcg_content(reference_id, start, end)
		values=list(merged_atcgdict.values())
		return list(zip(list(merged_atcgdict.keys()), list(map(lambda x:x[0], values)), list(map(lambda x:x[1], values))))

	def run_editing_percentage(self, call_editing_sites=False):
		result_dict = dict.fromkeys(self._chrom_sizes.keys(),list())
		for chrom,length in self._chrom_sizes.items():
			print(GET_CUR_TIME("Start analysing " + GET_BLUE("chromosome {}".format(chrom))))
			for i in tqdm.trange(0,length,self.frame_size):
				result_dict[chrom] = result_dict[chrom] + self.compute_nucleotides_coverage(chrom, i, i+self.frame_size)
				break
			break
		return result_dict

	def run(self):
		result_dict = dict.fromkeys(self._chrom_sizes.keys(),list())
		for chrom,length in self._chrom_sizes.items():
			print(GET_CUR_TIME("Start analysing " + GET_BLUE("chromosome {}".format(chrom))))
			for i in tqdm.trange(0,length,self.frame_size):
				result_dict[chrom] = result_dict[chrom] + self.call_editing_region(chrom, i, i+self.frame_size)
				break
			break
		return result_dict

	def run_par(self):
		result_dict = dict.fromkeys(self._chrom_sizes.keys(),list())
		for chrom,length in self._chrom_sizes.items():
			print(GET_CUR_TIME("Start analysing " + GET_BLUE("chromosome {}".format(chrom))))
			for i in tqdm.trange(0,length,self.frame_size):
				result_dict[chrom] = result_dict[chrom] + self.call_editing_region_par(chrom, i, i+self.frame_size)
		return result_dict

	def __call__(self):
		return self.run()

	