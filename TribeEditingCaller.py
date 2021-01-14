from .TribeReads import *
from .utils import *
from asyncio import Future, AbstractEventLoop, get_event_loop
from concurrent.futures import ThreadPoolExecutor
import asyncio
import threading
import numpy as np 
import scipy.stats as stats
from typing import Optional, Union, Tuple
import tqdm

class TribeCaller(object):
	def __init__(self, target_reads_path:str, 
				 control_reads_path:str, 
				 frame_size:int = 600000, 
				 bin_size:int = 100,
				 event_loop: Optional[AbstractEventLoop] = None,
				 num_threads:int=1,
				 run_async: bool = False):
		super(TribeCaller, self).__init__()
		self.target = TribeReads(target_reads_path)
		self.control = TribeReads(control_reads_path)
		if len(self.target._chrom_sizes) != len(self.control._chrom_sizes):
			print(get_cur_time("The length of chromosome number is not the same between target and control"))
		for i in self.target._chrom_sizes.keys():
			if self.target._chrom_sizes[i] != self.control._chrom_sizes[i]:
				print(get_cur_time("Different chromosome size found. Please check if use the same aligner"))
		self._chrom_sizes = self.target._chrom_sizes
		self.frame_size = frame_size
		self.bin_size = 100
		self.run_async_ = run_async
		self.event_loop: Optional[asyncio.AbstractEventLoop] = event_loop if event_loop else asyncio.get_event_loop()
		self.num_threads = num_threads

	def compute_ag_content(self, reference_id, start:int, end:int):
		target_agdict =  self.target.build_ag_dict(reference_id, start, end, self.bin_size)
		control_agdict = self.control.build_ag_dict(reference_id, start, end, self.bin_size)
		return target_agdict.merge(control_agdict, lambda x,y:(x,y))

	def call_editing_region(self,reference_id, start:int, end:int, threshold=3):
		merged_agdict = self.compute_ag_content(reference_id, start, end)
		def _frac(a,b):
			o = a/b if b != 0 else 0
			return o if o > 0 else 0
		log2exp = np.log2(list(map(lambda x:_frac(_frac(x[0][1],x[0][0]),_frac(x[1][1],x[1][0]))+1,merged_agdict.values())))
		# pval = list(map(lambda x: stats.hypergeom.sf(x[0][1], x[1][1]+x[1][0],x[0][1]+x[0][0],x[1][1]), merged_agdict.values()))
		return [list(merged_agdict.keys())[x] for x in range(len(log2exp)) if log2exp[x] > threshold]

	async def run_async(self, executor):
		future_dict = dict.fromkeys(self._chrom_sizes.keys(),list())
		result_dict = dict.fromkeys(self._chrom_sizes.keys(),list())
		for chrom,length in self._chrom_sizes.items():
			print(GET_CUR_TIME("Start analysing chromosome {}".format(chrom)))
			for i in tqdm.trange(0,length,self.frame_size):
				future_dict[chrom].append(self.event_loop.run_in_executor(executor, self.call_editing_region, chrom, i, i+self.frame_size))
				completed, pending = await asyncio.wait(future_dict[chrom])
				result_dict[chrom] = FLATTEN([t.result() for t in completed])
		return result_dict

	def __call__(self):
		with ThreadPoolExecutor(max_workers=self.num_threads) as executor:
			return self.run_async(executor)
