# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#

class TribeCriteria:
	def __init__(self, odds_threshold=0.5,
				       reads_coverage_threshold = 10, 
				       content_threshold=0.9, 
				       diff_threshold = 0.005,
				       pvalue_cutoff=0.05,
				       paired=False,
				       exclude_gap=False):
		self.odds_threshold = odds_threshold
		self.reads_coverage_threshold = reads_coverage_threshold
		self.content_threshold = content_threshold
		self.diff_threshold = diff_threshold
		self.pvalue_cutoff = pvalue_cutoff
		self.paired = paired
		self.exclude_gap = exclude_gap

	def render_args(self):
		return self.__dict__

	def set_args(self,arg_name,arg_value):
		return self.__setattr__(arg_name,arg_value)
