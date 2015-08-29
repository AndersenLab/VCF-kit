from cyvcf2 import VCF as cyvcf2
from collections import OrderedDict, deque
from itertools import islice
import re
import numpy as np


class vcf(cyvcf2):
    def __init__(self, filename):
        cyvcf2.__init__(self, filename)
        self.filename = filename

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

    def window(self, window_size, shift_method):
        """
            Generates windows of VCF data by positions or SNP
        """
        result_list = variant_interval(window_size = window_size, shift_method = shift_method)
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
                    while True:
                        result_list.append(line)
                        positions = result_list.positions()
                        max_pos = max(positions)
                        min_pos = min(positions)
                        while (min_pos >= result_list.upper_bound):
                                result_list.iterate_interval()
                                positions = result_list.positions()
                                min_pos = min(positions)
                        if not result_list.unique_chroms():
                            # Chromosome Reset
                            max_pos = 0
                            yield result_list[0:len(result_list)-1]
                            result_list.get_last()
                        elif max_pos >= result_list.lower_bound and max_pos > result_list.upper_bound and line not in result_list:
                            # If in current interval, do nothing.
                            pass
                        elif max_pos < result_list.lower_bound:
                            # If beneath interval - reset to beginning.
                            result_list.lower_bound = 0
                            result_list.upper_bound = window_size
                        elif max_pos >= result_list.upper_bound:
                            # If past interval, iterate and yield result.
                            if len(result_list) > 0:
                                yield result_list[0:len(result_list)-1]
                                result_list.get_last()
                                result_list.iterate_interval()
                            else:
                                pass

                                
                        line = self.next()

                elif shift_method == "POS-sliding":
                    if len(result_list) == 0:
                        result_list.append(line)
                    while True:
                        positions = [x.POS for x in result_list]
                        max_pos = max(positions)
                        min_pos = min(positions)
                        chrom_num = len(set([o.CHROM for o in result_list]))
                        if (max_pos - min_pos) > window_size or chrom_num > 1:
                            result_list.popleft()
                        else:
                            break
                    yield result_list
                    result_list.append(line)

        except StopIteration:
            if len(result_list) > 0:
                yield result_list


class variant_interval(deque):
    def __init__(self, interval = [], window_size = None, shift_method = None, varlist=None, lower_bound = 0):
        self.shift_method = shift_method
        self.window_size = window_size
        self.lower_bound = lower_bound
        self.upper_bound = self.lower_bound + window_size
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
        else:
            #CHROM = self[0].CHROM
            return ("", self.lower_bound, self.upper_bound)

    def iterate_interval(self):
        self.lower_bound += self.window_size
        self.upper_bound += self.window_size
        return self

    def __getitem__(self, index):
        if isinstance(index, slice):
            cut = islice(self, index.start, index.stop, index.step)
            return variant_interval(interval = cut, 
                                 window_size = self.window_size,
                                 shift_method = self.shift_method,
                                 lower_bound = self.lower_bound)
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


class variant_set:
    def __init__(self, varset):
        # self.shift_method = shift_method
        self.gt = np.vstack([x.gt_types for x in varset])
        self.segregating_sites = sum([np.any(p) for p in np.equal(3, x)])
        print self.segregating_sites


x = vcf("../test.vcf.gz")

for i in x.window(window_size=100000, shift_method="POS-Interval"):
    print(i)

    # var = variant_set(i)
    # break
