from .TribeReads import *
from .utils import *
from asyncio import Future, AbstractEventLoop, get_event_loop
import asyncio
import threading
import numpy as np 
import scipy.stats as stats
from typing import Optional, Union, Tuple
import tqdm
from concurrent.futures import ProcessPoolExecutor   

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
			print(get_cur_time("The length of chromosome number is not the same between target and control"))
		for i in self.target._chrom_sizes.keys():
			if self.target._chrom_sizes[i] != self.control._chrom_sizes[i]:
				print(get_cur_time("Different chromosome size found. Please check if use the same aligner"))
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

	def call_editing_region(self,reference_id, start:int, end:int, threshold=2):
		merged_asdict = self.compute_as_content(reference_id, start, end)
		def _frac(a,b):
			o = a/b if b != 0 else 0
			return o if o > 0 else 0
		res = list(map(lambda x: stats.fisher_exact([[x[0][1], x[1][1]],[x[0][0], x[1][0]]],"greater"), merged_asdict.values()))
		odds, pval = np.array(list(map(lambda x:x[0], res))),np.array(list(map(lambda x:x[1], res)))
		return [(list(merged_asdict.keys())[x],odds[x],pval[x]) for x in range(len(odds)) if odds[x] > threshold and pval[x] < 0.05]

	def run(self):
		result_dict = dict.fromkeys(self._chrom_sizes.keys(),list())
		for chrom,length in self._chrom_sizes.items():
			print(GET_CUR_TIME("Start analysing chromosome {}".format(chrom)))
			for i in tqdm.trange(0,length,self.frame_size):
				result_dict[chrom] = self.call_editing_region(chrom, i, i+self.frame_size)
		return result_dict

	def __call__(self):
		with ThreadPoolExecutor(max_workers=self.num_threads) as executor:
			return self.run_async(executor)
