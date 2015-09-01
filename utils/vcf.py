from cyvcf2 import VCF as cyvcf2
from collections import OrderedDict, deque
from itertools import islice, combinations
import re
import numpy as np
from copy import deepcopy
np.set_printoptions(threshold=np.nan)

class vcf(cyvcf2):
    def __init__(self, filename):
        cyvcf2.__init__(self, filename)
        self.filename = filename
        self.n = len(self.samples) # Number of Samples

        # Meta Data
        comp = re.compile(r'''^##(?P<key>[^<#]+?)=(?P<val>[^<#]+$)''', re.M)
        self.metadata = OrderedDict(comp.findall(self.raw_header))

        # Contigs
        self.contigs = dict(zip(
            re.compile("##contig=<ID=(.*?),").findall(self.raw_header),
            map(int, re.compile("##contig.*length=(.*?)>").findall(self.raw_header))
        ))
        # Info
        r = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.info_set = {x["id"]: x for x in [m.groupdict() for m in r.finditer(self.raw_header)]}

        # Filter
        r = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.filter_set = {x["id"]: x for x in [m.groupdict() for m in r.finditer(self.raw_header)]}


        # Format
        r = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
        self.format_set = {x["id"]: x for x in [m.groupdict() for m in r.finditer(self.raw_header)]}

    def window(self, shift_method, window_size, step_size = None):
        """
            Generates windows of VCF data by positions or SNP
        """
        result_list = variant_interval(window_size = window_size, step_size = step_size, shift_method = shift_method)
        try:
            while True:
                line = self.next()
                #
                # SNP-Sliding and SNP-Interval
                #
                if shift_method in ["SNP-Sliding", "SNP-Interval"]:
                    result_list.append(line)
                    if not result_list.unique_chroms():
                        yield result_list[0:len(result_list)-1]
                        result_list.get_last()
                    elif len(result_list) == window_size:
                        yield result_list
                        if shift_method == "SNP-Interval":
                            result_list.clear()
                #
                # POS-Interval
                #
                elif shift_method == "POS-Interval":
                    result_list.append(line)
                    last_chrom = line.CHROM
                    chrom = line.CHROM
                    while True:
                        if chrom != last_chrom:
                            yield result_list.filter_within_bounds()
                            result_list.get_last()
                            result_list.lower_bound = 0
                            result_list.upper_bound = window_size
                            chrom = line.CHROM
                        if result_list.lower_bound  <= line.POS < result_list.upper_bound:
                            line = self.next()
                            result_list.append(line)
                        if line.POS < result_list.lower_bound:
                            result_list.iterate_interval()
                        elif line.POS >= result_list.upper_bound:
                            yield result_list.filter_within_bounds()
                            result_list.iterate_interval()
                        last_chrom = line.CHROM
                        #print line.POS
                #
                # POS-Sliding
                #
                elif shift_method == "POS-Sliding":
                    result_list.append(line)
                    while True:
                        positions = result_list.positions()
                        max_pos = max(positions)
                        min_pos = min(positions)
                        chrom_num = len(set([o.CHROM for o in result_list]))
                        if (max_pos - min_pos) > window_size or chrom_num > 1:
                            result_list.popleft()
                        else:
                            break
                    yield result_list

        except StopIteration:
            if len(result_list) > 0:
                if shift_method == "POS-Interval":
                    yield result_list.filter_within_bounds()
                else:
                    yield result_list

class variant_interval(deque):
    def __init__(self, interval = [], window_size = None, step_size = None, shift_method = None, varlist=None, lower_bound = 0):
        self.shift_method = shift_method
        self.window_size = window_size
        self.step_size = step_size
        if self.step_size == None:
            self.step_size = window_size
        self.lower_bound = lower_bound
        self.upper_bound = self.lower_bound + window_size - 1
        if shift_method in ["SNP-Sliding", "SNP-Interval"]:
            super(variant_interval, self).__init__(interval, maxlen = self.window_size)
        else:
            super(variant_interval, self).__init__(interval)

    def positions(self):
        return [x.POS for x in self]

    def interval(self):
        if self.shift_method in ["SNP-Sliding", "SNP-Interval"]:
            POS = self.positions()
            CHROM = self[0].CHROM
            return (CHROM, min(POS), max(POS))
        elif self.shift_method == "POS-Sliding":
            # Space interval evenly around items
            mean_pos = np.mean(self.positions())
            half_window = self.window_size/2
            lower_bound = int(mean_pos - half_window)
            CHROM = self[0].CHROM
            if lower_bound < 0:
                lower_bound = 0
            upper_bound = int(mean_pos + half_window)
            return (CHROM, lower_bound, upper_bound)
        else:
            CHROM = self[0].CHROM
            return (CHROM, self.lower_bound, self.upper_bound)

    def iterate_interval(self):
        self.lower_bound += self.step_size
        self.upper_bound += self.step_size
        return self

    def filter_within_bounds(self):
        cut = [x for x in self if x.POS >= self.lower_bound and x.POS < self.upper_bound]
        return variant_interval(interval = cut, 
                             window_size = self.window_size,
                             shift_method = self.shift_method,
                             lower_bound = self.lower_bound,
                             step_size = self.step_size)


    def __getitem__(self, index):
        if isinstance(index, slice):
            cut = islice(self, index.start, index.stop, index.step)
            return variant_interval(interval = cut, 
                                 window_size = self.window_size,
                                 shift_method = self.shift_method,
                                 lower_bound = self.lower_bound,
                                 step_size = self.step_size)
        return deque.__getitem__(self, index)


    def get_last(self):
        """
            Removes all items from the variant interval and
            sets to last item
        """
        super(variant_interval, self).__init__(self[len(self)-1:len(self)], maxlen = self.window_size)
        

    def unique_chroms(self):
        """
            Checks to see whether the variants
            all belong to the same chromosome
        """
        return len(set([var.CHROM for var in self])) == 1 and len(self) > 0

    def __repr__(self):
        formatted_variants = ["{chrom}:{pos}".format(chrom=var.CHROM, pos=var.POS) for var in self]
        formatted_variants = "\t".join(formatted_variants)
        interval = self.interval()
        if interval is not None:
            return "{0}:{1}-{2}".format(*self.interval()) + " -> " + formatted_variants
        else:
            return  str(formatted_variants) + " Format"

