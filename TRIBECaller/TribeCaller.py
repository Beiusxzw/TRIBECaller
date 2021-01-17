# -*- coding: utf-8 -*-
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

from asyncio import Future, AbstractEventLoop, get_event_loop
import asyncio
import threading
import numpy as np 
import scipy.stats as stats
from typing import Optional, Union, Tuple
import tqdm
from concurrent.futures import ProcessPoolExecutor   

# ------------------------- #
# Main Class
# ------------------------- #

class TribeCaller(object):
	def __init__(self, target_reads_path:str, 
				 control_reads_path:str, 
				 frame_size:int = 600000, 
				 bin_size:int = 1,
				 event_loop: Optional[AbstractEventLoop] = None,
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

		self.run_async_ = run_async
		self.event_loop: Optional[asyncio.AbstractEventLoop] = event_loop if event_loop else asyncio.get_event_loop()
		self.num_threads = num_threads

	def compute_as_content(self, reference_id, start:int, end:int):
		target_asdict =  self.target.build_as_dict(reference_id, start, end, self.bin_size)
		control_asdict = self.control.build_as_dict(reference_id, start, end, self.bin_size)
		return target_asdict.merge(control_asdict, lambda x,y:(x,y))

	def compute_region_as_percentage(self,reference_id,start:int, end:int=None,call_editing_sites=False):
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
		result = []
		for i in range(start,end):
			t,c = target_atcgdict[i],control_atcgdict[i]
			result.append([[t[0]/t[4], t[1]/t[4],t[2]/t[4],t[3]/t[4]] if t[4] > 0 else [0,0,0,0] ,[c[0]/c[4], c[1]/c[4],c[2]/c[4],c[3]/c[4]] if c[4] > 0 else [0,0,0,0]])
		if call_editing_sites:
			res = list(map(lambda x: stats.fisher_exact([[x[0][2]+x[0][3], x[1][2]+x[1][3]],[x[0][0], x[1][0]]], "greater"), MERGE_DICT(target_atcgdict, control_atcgdict, lambda x,y:(x,y)).values()))
			odds, pval = np.array(list(map(lambda x:x[0], res))),np.array(list(map(lambda x:x[1], res)))
			return list(range(start,end)),result,odds,pval
		return list(range(start,end)),result

	def call_editing_region(self,reference_id, start:int, end:int, threshold=2, test="fisher-exact"):
		merged_asdict = self.compute_as_content(reference_id, start, end)
		values=list(merged_asdict.values())
		res = list(map(lambda x: stats.fisher_exact([[x[0][1], x[1][1]],[x[0][0], x[1][0]]], "greater"), values))
		odds, pval = np.array(list(map(lambda x:x[0], res))),np.array(list(map(lambda x:x[1], res)))
		return [(list(merged_asdict.keys())[x],odds[x],pval[x]) for x in range(len(odds)) if odds[x] > threshold and pval[x] < 0.05 and values[x][1][0]/values[x][1][2] > 0.8]

	def call_editing_region_coverage(self,reference_id, start:int, end:int, threshold=2, test="fisher-exact"):
		print(GET_CUR_TIME("Computing bam coverage and nucleotides content"))
		if end - start > 100000:
			raise TypeError("the region length should not exceed 100000bp")
		merged_asdict = self.compute_as_content(reference_id, start, end)
		values=list(merged_asdict.values())
		res = list(map(lambda x: stats.fisher_exact([[x[0][1], x[1][1]],[x[0][0], x[1][0]]], "greater"), values))
		odds, pval = np.array(list(map(lambda x:x[0], res))),np.array(list(map(lambda x:x[1], res)))
		with np.errstate(invalid='ignore'):
			return list(merged_asdict.keys()), list(map(lambda y: (y[0][2],y[1][2]), merged_asdict.values())), (odds > threshold) & (pval < 0.05) & (np.array(list(map(lambda x:x[1][0]/x[1][2] > 0.8,values))))

	def run(self):
		result_dict = dict.fromkeys(self._chrom_sizes.keys(),list())
		for chrom,length in self._chrom_sizes.items():
			print(GET_CUR_TIME("Start analysing " + GET_BLUE("chromosome {}".format(chrom))))
			for i in tqdm.trange(0,length,self.frame_size):
				result_dict[chrom] = result_dict[chrom] + self.call_editing_region(chrom, i, i+self.frame_size)
		return result_dict

	def __call__(self):
		return self.run()

	