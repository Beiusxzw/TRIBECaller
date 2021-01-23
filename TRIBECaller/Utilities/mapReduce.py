# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
# Description: Functions and Methods for functional programming
#
# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#

#from pathos.multiprocessing import Pool
from multiprocessing import Process, Manager, Queue
import time

class MapReduce(object):
	def __init__(self,  map_func, reduce_func, chunk_size, num_workers=1):
		super(MapReduce, self).__init__()
		self.map_func = map_func
		self.reduce_func = reduce_func
		self.chunk_size = chunk_size
		self.Pool = Pool(num_workers)

	def partition(self, iterable):
		print("total length: {}".format(len(iterable)))
		return [iterable[x:x+self.chunk_size] for x in range(0, len(iterable), self.chunk_size)]

	def sync_map_reduce(self, iterable):
		"""
		For testing map_func and reduce_func
		@args iterable to be mapped
		@return reduced result, synchronously
		"""
		if len(iterable) > self.chunk_size:
			return self.reduce_func(self.sync_map_reduce(iterable[:int(len(iterable)/2)]),
									self.sync_map_reduce(iterable[int(len(iterable)/2):]))
		else:
			return self.map_func(iterable)

	def par_map_reduce(self, iterable):
		"""
		parallel implementation of mapReduce using map method from multiprocessing.Pool
		@args iterable to be mapped
		@returns the reduce result of the input iterable
		"""
		iterable = self.partition(iterable)
		print("partition: {}".format(len(iterable)))
		start = time.time()
		result = self.Pool.map(self.map_func, iterable)
		print("Pool maped: {}".format(time.time()-start))
		return self.reduce_func(result)

	def map_func_wrap(self, procnum, return_dict, *args):
		return_dict[procnum] = self.map_func(*args)

	def proc_map_reduce(self, iterable):
		iterable = list(self.partition(iterable))
		print("partition: {}".format(len(iterable)))
		manager = Manager()
		return_dict = manager.dict()
		proc_queue = []
		for proc,i in enumerate(iterable):
			p = Process(target=self.map_func_wrap, args=(proc, return_dict,i))
			p.Daemon = True
			p.start()
			proc_queue.append(p)
		for p in proc_queue:
			p.join()
		while True:
			if any(proc.is_alive() for proc in proc_queue):
				time.sleep(1)
			else:
				return self.reduce_func(return_dict.values())