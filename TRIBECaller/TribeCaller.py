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
import gzip
from functools import partial

# ------------------------- #
# Main Class
# ------------------------- #


class TribeCaller(object):
    def __init__(self, target_reads_path: str,
                 control_reads_path: str,
                 frame_size: int = 500000,
                 bin_size: int = 1,
                 num_threads: int = 1):
        """
        class TribeCaller: main functional class for TRIBECaller.
        @args target_reads_path: str, path to the target bam file
        @args control reads path: str, path to the control bam file
        @args frame size: Framesize for computing the atcg dictionary
        @args bin_size: int, resolution of the editing sites. default=1
        @args num_threads:int, number of threads to use in the main function
        """
        super(TribeCaller, self).__init__()
        self.target = TribeReads(target_reads_path)
        self.control = TribeReads(control_reads_path)
        if len(self.target._chrom_sizes) != len(self.control._chrom_sizes):
            print(GET_CUR_TIME(
                "The length of chromosome number is not the same between target and control"))
        for i in self.target._chrom_sizes.keys():
            if self.target._chrom_sizes[i] != self.control._chrom_sizes[i]:
                print(GET_CUR_TIME(
                    "Different chromosome size found. Please check if use the same aligner"))
        self._chrom_sizes = self.target._chrom_sizes
        self.frame_size = frame_size
        self.bin_size = bin_size
        self.num_threads = num_threads

    def compute_atcg_content(self, reference_id, start: int, end: int, paired=False, exclude_gap=False):
        """
        compute A,T,C,G content in each nucleotides position
        @args reference_id: contig name of interest
        @args start: integer, start position
        @args end: integer, end position
        """
        target_atcgdict = self.target.build_atcg_dict(
            reference_id, start, end, self.bin_size, paired=paired, exclude_gap=exclude_gap)
        control_atcgdict = self.control.build_atcg_dict(
            reference_id, start, end, self.bin_size, paired=paired, exclude_gap=exclude_gap)
        return target_atcgdict.merge(control_atcgdict, lambda x, y: (x, y))

    def compute_atcg_content_par(self, reference_id, start_end_list, n_threads, paired=False, exclude_gap=False):
        """
        parallel implementation of compute_atcg_content
        @args reference_id: contig name of interest
        @args start_end_list: list consist of start-end positions
        @args n_threads: integer indicating number of threads to use. 
        """
        target_atcgdict = self.target.build_atcg_dict_par(
            reference_id, start_end_list, n_threads, self.bin_size, paired=paired, exclude_gap=exclude_gap)
        result = []
        control_atcgdict = self.control.build_atcg_dict_par(
            reference_id, start_end_list, n_threads, self.bin_size, paired=paired, exclude_gap=exclude_gap)
        for i in zip(target_atcgdict, control_atcgdict):
            result.append(i[0].merge(i[1], lambda x, y: (x, y)))
        return result

    def compute_region_as_percentage(self, reference_id,
                                     start: int,
                                     end: int,
                                     call_editing_sites,
                                     odds_threshold,
                                     content_threshold,
                                     reads_coverage_threshold,
                                     pvalue_cutoff,
                                     diff_threshold,
                                     paired,
                                     exclude_gap):
        """
        for use of plotEditingRegion.py
        compute A,T,C,G content in each nucleotides position
        @args reference_id: contig name of interest
        @args start: integer, start position
        @args end: integer, end position
        @args odds_threshold: 
        @args content_threshold
        @args reads_coverage_threshold
        @args pvalue_cutoff
        @args diff_threshold
        @args paired
        @args exclude gap
        """
        def _helper(x):
            return 0 if x < 0 else x
        if end:
            if end - start > 200:
                raise TypeError("the region length should not exceed 200bp")
        else:
            start = start - 20
            end = start + 40
        target_atcgdict = self.target.build_atcg_dict(
            reference_id, start, end, self.bin_size, paired=paired, exclude_gap=exclude_gap)
        target_atcgdict = {k: target_atcgdict[k] if k in target_atcgdict.keys() else (
            0, 0, 0, 0, 0) for k in range(start, end)}
        control_atcgdict = self.control.build_atcg_dict(
            reference_id, start, end, self.bin_size, paired=paired, exclude_gap=exclude_gap)
        control_atcgdict = {k: control_atcgdict[k] if k in control_atcgdict.keys(
        ) else (0, 0, 0, 0, 0) for k in range(start, end)}
        count, result = [], []
        for i in range(start, end):
            t, c = target_atcgdict[i], control_atcgdict[i]
            if call_editing_sites:
                result.append([[t[0]/t[4], t[1]/t[4], t[2]/t[4], t[3]/t[4], t[4]] if t[4] > 0 else [0, 0, 0, 0, 0], [
                              c[0]/c[4], c[1]/c[4], c[2]/c[4], c[3]/c[4], c[4]] if c[4] > 0 else [0, 0, 0, 0, 0]])
                count.append([[t[0], t[1], t[2], t[3], t[4]] if t[4] > 0 else [0, 0, 0, 0, 0], [
                             c[0], c[1], c[2], c[3], c[4]] if c[4] > 0 else [0, 0, 0, 0, 0]])
            else:
                result.append([[t[0]/t[4], t[1]/t[4], t[2]/t[4], t[3]/t[4]] if t[4] > 0 else [0, 0, 0, 0], [
                              c[0]/c[4], c[1]/c[4], c[2]/c[4], c[3]/c[4]] if c[4] > 0 else [0, 0, 0, 0]])
        if call_editing_sites:
            diff = list(map(lambda x: (_helper((x[0][3]/(x[0][0]+x[0][3])) - (x[1][3]/(x[1][0]+x[1][3]))) if (x[0][0]+x[0][3] > 0 and x[1][0]+x[1][3] > 0) else 0, _helper(
                (x[0][2]/(x[0][1]+x[0][2]) - (x[1][2]/(x[1][1]+x[1][2])))) if (x[0][1]+x[0][2] > 0 and x[1][1]+x[1][2] > 0) else 0) if x[0][4] > 5 and x[1][4] > 5 else 0, count))
            res = list(map(lambda x: (stats.fisher_exact([[x[0][3], x[1][3]], [x[0][0], x[1][0]]], "greater"), stats.fisher_exact(
                [[x[0][2], x[1][2]], [x[0][1], x[1][1]]], "greater")), count))
            odds_f, pval_f, odds_r, pval_r = np.array(list(map(lambda x: x[0], map(lambda x: x[0], res)))), np.array(list(map(lambda x: x[1], map(
                lambda x: x[0], res)))), np.array(list(map(lambda x: x[0], map(lambda x: x[1], res)))), np.array(list(map(lambda x: x[1], map(lambda x: x[1], res))))
            with np.errstate(invalid='ignore'):
                return list(range(start, end)), list(map(lambda x: [x[0][:4], x[1][:4]], result)), ((odds_f > odds_threshold) & (pval_f < pvalue_cutoff) & (np.array(list(map(lambda x: x[1][0]/x[1][4] >= content_threshold if x[1][4] > reads_coverage_threshold else False, count)))) & (np.array(list(map(lambda x: x[0] if type(x) == tuple else x, diff))) > diff_threshold) ) | ((odds_r > odds_threshold) & (pval_r < pvalue_cutoff) & (np.array(list(map(lambda x: x[1][1]/x[1][4] >= content_threshold if x[1][4] > reads_coverage_threshold else False, count))))) & (np.array(list(map(lambda x: x[1] if type(x) == tuple else x, diff))) > diff_threshold)
        return list(range(start, end)), result

    def call_editing_region_coverage(self, reference_id,
                                     start: int,
                                     end: int,
                                     odds_threshold,
                                     content_threshold,
                                     reads_coverage_threshold,
                                     pvalue_cutoff,
                                     diff_threshold,
                                     paired,
                                     exclude_gap,
                                     output_complex=True):
        """
        for use of plotEditingRegion.py
        compute read coverage of every nucleotide poition, and call editing sites
        @args reference_id: contig name of interest
        @args start: integer, start position
        @args end: integer, end position
        @args odds_threshold: 
        @args content_threshold
        @args reads_coverage_threshold
        @args pvalue_cutoff
        @args diff_threshold
        @args paired
        @args exclude gap
        """
        def _helper(x):
            return 0 if x < 0 else x
        print(GET_CUR_TIME("Computing bam coverage and nucleotides content"))
        # if end - start > 200000:
        #	raise TypeError("the region length should not exceed 200000bp")
        merged_atcgdict = self.compute_atcg_content(
            reference_id, start, end, paired=paired, exclude_gap=exclude_gap)
        values = list(merged_atcgdict.values())
        diff = list(map(lambda x: (_helper((x[0][3]/(x[0][0]+x[0][3])) - (x[1][3]/(x[1][0]+x[1][3]))) if (x[0][0]+x[0][3] > 0 and x[1][0]+x[1][3] > 0) else 0, _helper(
            (x[0][2]/(x[0][1]+x[0][2]) - (x[1][2]/(x[1][1]+x[1][2])))) if (x[0][1]+x[0][2] > 0 and x[1][1]+x[1][2] > 0) else 0) if x[0][4] > 5 and x[1][4] > 5 else 0, values))
        res = list(map(lambda x: (stats.fisher_exact([[x[0][3], x[1][3]], [x[0][0], x[1][0]]], "greater"), stats.fisher_exact(
            [[x[0][2], x[1][2]], [x[0][1], x[1][1]]], "greater")), values))
        odds_f, pval_f, odds_r, pval_r = np.array(list(map(lambda x: x[0], map(lambda x: x[0], res)))), np.array(list(map(lambda x: x[1], map(
            lambda x: x[0], res)))), np.array(list(map(lambda x: x[0], map(lambda x: x[1], res)))), np.array(list(map(lambda x: x[1], map(lambda x: x[1], res))))
        with np.errstate(invalid='ignore'):
            return list(merged_atcgdict.keys()), list(map(lambda y: (y[0][4], y[1][4]), merged_atcgdict.values())), ((odds_f > odds_threshold) & (pval_f < pvalue_cutoff) & (np.array(list(map(lambda x: x[1][0]/x[1][4] >= content_threshold and x[0][4] > reads_coverage_threshold and x[1][4] > reads_coverage_threshold, values)))) & (np.array(list(map(lambda x: x[0] if type(x) == tuple else x, diff))) > diff_threshold)), ((odds_r > odds_threshold) & (pval_r < pvalue_cutoff) & (np.array(list(map(lambda x: x[1][1]/x[1][4] >= content_threshold and x[0][4] > reads_coverage_threshold and x[1][4] > reads_coverage_threshold, values)))) & (np.array(list(map(lambda x: x[1] if type(x) == tuple else x, diff))) > diff_threshold)), list(map(lambda x: max(x) if type(x) == tuple else x, diff))

    def call_editing_region(self, reference_id,
                            start: int,
                            end: int,
                            odds_threshold,
                            content_threshold,
                            reads_coverage_threshold,
                            pvalue_cutoff,
                            diff_threshold,
                            paired,
                            exclude_gap,
                            output_complex=True):
        """
        call editing events in a certain region
        @args reference_id: contig name of interest
        @args start: integer, start position
        @args end: integer, end position
        @args odds_threshold: 
        @args content_threshold
        @args reads_coverage_threshold
        @args pvalue_cutoff
        @args diff_threshold
        @args paired
        @args exclude gap
        """
        merged_atcgdict = self.compute_atcg_content(
            reference_id, start, end, paired=paired, exclude_gap=exclude_gap)
        return self.call_editing_region_wrap(merged_atcgdict, odds_threshold, content_threshold, reads_coverage_threshold, pvalue_cutoff, diff_threshold, paired, exclude_gap, output_complex)

    def call_editing_region_wrap(self, merged_atcgdict,
                                 odds_threshold,
                                 content_threshold,
                                 reads_coverage_threshold,
                                 pvalue_cutoff,
                                 diff_threshold,
                                 paired,
                                 exclude_gap,
                                 output_complex=True):
        """
        parallel implementation of calling editing events in a certain region
        @args reference_id: contig name of interest
        @args start: integer, start position
        @args end: integer, end position
        @args odds_threshold: 
        @args content_threshold
        @args reads_coverage_threshold
        @args pvalue_cutoff
        @args diff_threshold
        @args paired
        @args exclude gap
        """
        values = list(merged_atcgdict.values())
        diff = list(map(lambda x: ((x[0][3]/(x[0][0]+x[0][3])) - (x[1][3]/(x[1][0]+x[1][3])) if (x[0][0]+x[0][3] > 0 and x[1][0]+x[1][3] > 0)
                                   else 0, (x[0][2]/(x[0][1]+x[0][2]) - (x[1][2]/(x[1][1]+x[1][2]))) if (x[0][1]+x[0][2] > 0 and x[1][1]+x[1][2] > 0) else 0), values))
        res = list(map(lambda x: (stats.fisher_exact([[x[0][3], x[1][3]], [x[0][0], x[1][0]]], "greater"), stats.fisher_exact(
            [[x[0][2], x[1][2]], [x[0][1], x[1][1]]], "greater")), values))
        odds_f, pval_f, odds_r, pval_r = np.array(list(map(lambda x: x[0], map(lambda x: x[0], res)))), np.array(list(map(lambda x: x[1], map(
            lambda x: x[0], res)))), np.array(list(map(lambda x: x[0], map(lambda x: x[1], res)))), np.array(list(map(lambda x: x[1], map(lambda x: x[1], res))))
        if output_complex:
            out = []
            for x in range(len(res)):
                if (odds_f[x] > odds_threshold and pval_f[x] < pvalue_cutoff and values[x][1][0]/values[x][1][4] >= content_threshold and values[x][0][4] > reads_coverage_threshold and values[x][1][4] > reads_coverage_threshold and diff[x][0] > diff_threshold):
                    out.append((list(merged_atcgdict.keys())[
                               x], diff[x][0], values[x][0][0], values[x][0][3], values[x][1][0], values[x][1][3], pval_f[x], '+'))
                if (odds_r[x] > odds_threshold and pval_r[x] < pvalue_cutoff and values[x][1][1]/values[x][1][4] >= content_threshold and values[x][0][4] > reads_coverage_threshold and values[x][1][4] > reads_coverage_threshold and diff[x][1] > diff_threshold):
                    out.append((list(merged_atcgdict.keys())[
                               x], diff[x][1], values[x][0][1], values[x][0][2], values[x][1][1], values[x][1][2], pval_r[x], '-'))
            return out
        else:
            out = []
            for x in range(len(res)):
                if (odds_f[x] > odds_threshold and pval_f[x] < pvalue_cutoff and values[x][1][0]/values[x][1][4] >= content_threshold and values[x][0][4] > reads_coverage_threshold and values[x][1][4] > reads_coverage_threshold and diff[x][0] > diff_threshold):
                    out.append((list(merged_atcgdict.keys())[
                               x], diff[x][0], pval_f[x], '+'))
                if (odds_r[x] > odds_threshold and pval_r[x] < pvalue_cutoff and values[x][1][1]/values[x][1][4] >= content_threshold and values[x][0][4] > reads_coverage_threshold and values[x][1][4] > reads_coverage_threshold and diff[x][1] > diff_threshold):
                    out.append((list(merged_atcgdict.keys())[
                               x], diff[x][1], pval_r[x], '-'))
            return out

    def compute_nucleotides_coverage(self, reference_id, start: int, end: int, paired=False, exclude_gap=False):
        """
        compute nucleotides coverage in a certain region. Do not call editing events.
        """
        merged_atcgdict = self.compute_atcg_content(
            reference_id, start, end, paired=paired, exclude_gap=exclude_gap)
        return self.compute_nucleotides_coverage_wrap(merged_atcgdict)

    def compute_nucleotides_coverage_wrap(self, merged_atcgdict):
        """
        Parallel implementation for computing nucleotides coverage in a certain region. Do not call editing events.
        """
        values = list(merged_atcgdict.values())
        return list(zip(list(merged_atcgdict.keys()), list(map(lambda x: x[0], values)), list(map(lambda x: x[1], values))))

    def write_editing_events(self, fd, chrom, data):
        for i in data:
            for j in i:
                conc = list(map(str, j[1:]))
                fd.write(
                    "\t".join([chrom] + [str(j[0]), str(j[0]+1)] + conc) + "\n")

    def write_nucleotides_percentage(self, fd, chrom, data):
        for i in data:
            for j in i:
                conc = list(map(str, j[1:]))
                fd.write("\t".join([chrom] + [str(j[0]), str(j[0]+1)] +
                                   list(map(str, j[1])) + list(map(str, j[2]))) + "\n")

    def run(self, criteria, out_prefix, g_zip=False, contig=None):
        f = gzip.open(out_prefix + ".bed.gz",
                      "wt") if g_zip else open(out_prefix + ".txt", "w+")
        if contig:
            if type(contig) == str:
                contig = [contig]
        if contig:
            for chrom, length in self._chrom_sizes.items():
                if chrom in contig:
                    print(GET_CUR_TIME("Start analysing " +
                                       GET_BLUE("chromosome {}".format(chrom))))
                    for i in tqdm.trange(0, length, self.frame_size):
                        for j in self.call_editing_region(chrom, i, i+self.frame_size, **criteria.render_args()):
                            conc = list(map(str, j[1:]))
                            f.write(
                                "\t".join([chrom] + [str(j[0]), str(j[0]+1)] + conc) + "\n")
        else:
            for chrom, length in self._chrom_sizes.items():
                print(GET_CUR_TIME("Start analysing " +
                                   GET_BLUE("chromosome {}".format(chrom))))
                for i in tqdm.trange(0, length, self.frame_size):
                    for j in self.call_editing_region(chrom, i, i+self.frame_size, **criteria.render_args()):
                        conc = list(map(str, j[1:]))
                        f.write(
                            "\t".join([chrom] + [str(j[0]), str(j[0]+1)] + conc) + "\n")
        f.close()

    def run_par(self, criteria, out_prefix, n_threads, g_zip=False, contig=None):
        print(GET_CUR_TIME("Running program using " +
                           GET_BLUE("{} threads".format(n_threads))))
        map_func = partial(self.call_editing_region_wrap,
                           **criteria.render_args())
        map_reduce = MapReduce(
            map_func=map_func, reduce_func=None, num_workers=n_threads)
        f = gzip.open(out_prefix + ".bed.gz",
                      "wt") if g_zip else open(out_prefix + ".bed", "w+")
        if contig:
            if type(contig) == str:
                contig = [contig]
        temp_data = ThreadDataList(n_threads)
        if contig:
            for chrom, length in self._chrom_sizes.items():
                if chrom in contig:
                    map_reduce.reduce_func = partial(
                        self.write_editing_events, f, chrom)
                    print(GET_CUR_TIME("Start analysing " +
                                       GET_BLUE("chromosome {}".format(chrom))))
                    for i in tqdm.trange(0, length, self.frame_size*n_threads):
                        for j in range(i, i + self.frame_size * n_threads, self.frame_size):
                            if j + self.frame_size <= length:
                                temp_data.append((j, j+self.frame_size))
                        map_reduce.proc_map_reduce(self.compute_atcg_content_par(
                            chrom, temp_data, n_threads, paired=criteria.paired, exclude_gap=criteria.exclude_gap))
                        temp_data.clear()

        else:
            for chrom, length in self._chrom_sizes.items():
                map_reduce.reduce_func = partial(
                    self.write_editing_events, f, chrom)
                print(GET_CUR_TIME("Start analysing " +
                                   GET_BLUE("chromosome {}".format(chrom))))
                for i in tqdm.trange(0, length, self.frame_size*n_threads):
                    for j in range(i, i + self.frame_size * n_threads, self.frame_size):
                        if j + self.frame_size <= length:
                            temp_data.append((j, j+self.frame_size))
                    map_reduce.proc_map_reduce(self.compute_atcg_content_par(
                        chrom, temp_data, n_threads, paired=criteria.paired, exclude_gap=criteria.exclude_gap))
                    temp_data.clear()

    def run_editing_percentage(self, criteria, out_prefix, contig=None, g_zip=False):
        if contig:
            if type(contig) == str:
                contig = [contig]
        result_dict = dict.fromkeys(self._chrom_sizes.keys(), list())
        f = gzip.open(out_prefix + ".txt.gz",
                      "wt") if g_zip else open(out_prefix + ".txt", "w+")
        if contig:
            for chrom, length in self._chrom_sizes.items():
                if chrom in contig:
                    print(GET_CUR_TIME("Start analysing " +
                                       GET_BLUE("chromosome {}".format(chrom))))
                    for i in tqdm.trange(0, length, self.frame_size):
                        for j in self.compute_nucleotides_coverage(chrom, i, i+self.frame_size, paired=criteria.paired, exclude_gap=criteria.exclude_gap):
                            conc = list(map(str, j[1:]))
                            f.write("\t".join(
                                [chrom] + [str(j[0]), str(j[0]+1)] + list(map(str, j[1])) + list(map(str, j[2]))) + "\n")
        else:
            for chrom, length in self._chrom_sizes.items():
                print(GET_CUR_TIME("Start analysing " +
                                   GET_BLUE("chromosome {}".format(chrom))))
                for i in tqdm.trange(0, length, self.frame_size):
                    for j in self.compute_nucleotides_coverage(chrom, i, i+self.frame_size, paired=criteria.paired, exclude_gap=criteria.exclude_gap):
                        conc = list(map(str, j[1:]))
                        f.write("\t".join([chrom] + [str(j[0]), str(j[0]+1)] +
                                          list(map(str, j[1])) + list(map(str, j[2]))) + "\n")
        f.close()

    def run_editing_percentage_par(self, out_prefix,  n_threads, contig=None, g_zip=False, paired=False, exclude_gap=False):
        map_reduce = MapReduce(
            map_func=self.compute_nucleotides_coverage_wrap, reduce_func=None, num_workers=n_threads)
        temp_data = ThreadDataList(n_threads)
        if contig:
            if type(contig) == str:
                contig = [contig]
        result_dict = dict.fromkeys(self._chrom_sizes.keys(), list())
        f = gzip.open(out_prefix + ".txt.gz",
                      "wt") if g_zip else open(out_prefix + ".txt", "w+")
        if contig:
            for chrom, length in self._chrom_sizes.items():
                map_reduce.reduce_func = partial(
                    self.write_nucleotides_percentage, f, chrom)
                if chrom in contig:
                    print(GET_CUR_TIME("Start analysing " +
                                       GET_BLUE("chromosome {}".format(chrom))))
                    for i in tqdm.trange(0, length, self.frame_size * n_threads):
                        for j in range(i, i + self.frame_size * n_threads, self.frame_size):
                            if j + self.frame_size <= length:
                                temp_data.append((j, j+self.frame_size))
                        map_reduce.proc_map_reduce(self.compute_atcg_content_par(
                            chrom, temp_data, paired=paired, exclude_gap=exclude_gap))
                        temp_data.clear()

        else:
            for chrom, length in self._chrom_sizes.items():
                map_reduce.reduce_func = partial(
                    self.write_nucleotides_percentage, f, chrom)
                print(GET_CUR_TIME("Start analysing " +
                                   GET_BLUE("chromosome {}".format(chrom))))
                for i in tqdm.trange(0, length, self.frame_size*n_threads):
                    for j in range(i, i + self.frame_size * n_threads, self.frame_size):
                        if j + self.frame_size <= length:
                            temp_data.append((j, j+self.frame_size))
                    map_reduce.proc_map_reduce(self.compute_atcg_content_par(
                        chrom, temp_data, paired=paired, exclude_gap=exclude_gap))
                    temp_data.clear()
        f.close()
