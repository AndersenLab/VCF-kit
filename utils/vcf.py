from cyvcf2 import VCF as cyvcf2
from collections import OrderedDict, deque
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
                        map(int,re.compile("##contig.*length=(.*?)>").findall(self.raw_header))
                        ))
        # Info
        r = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.info_set = {x["id"]:x for x in [m.groupdict() for m in r.finditer(self.raw_header)]}
        
        # Filter
        r = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.filter_set = {x["id"]:x for x in [m.groupdict() for m in r.finditer(self.raw_header)]}
        

        # Format
        r = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
        self.format_set = {x["id"]:x for x in [m.groupdict() for m in r.finditer(self.raw_header)]}

    def window(self, windowsize, shift_method):
        """
            Generates windows of VCF data by positions or SNP
        """
        if shift_method == ["SNP","sliding"]:
            result_list = deque(maxlen = windowsize)
        else:
            result_list = deque()
        # Defaults
        curr_interval = [0, windowsize]
        try:
            while True:
                line = self.next()
                if shift_method == ["SNP","sliding"]:
                    result_list.append(line)
                    if len(set([o.CHROM for o in result_list])) > 1:
                        yield list(result_list)[:-1]
                        result_list = [result_list[-1]]
                    elif len(result_list) == windowsize:
                        yield result_list
                elif shift_method == ["SNP","interval"]:
                    result_list.append(line)
                    # Test for new chromosomes
                    if len(set([o.CHROM for o in result_list])) > 1:
                        yield result_list[:-1]
                        result_list = [result_list[-1]]
                    elif len(result_list) == windowsize:
                        yield result_list
                        result_list = []
                elif shift_method == ["POS", "interval"]:
                    result_list = []
                    while True:
                        result_list.append(line)
                        positions = [x.POS for x in result_list]
                        chrom_set = len(set([x.CHROM for x in result_list]))
                        max_pos = max(positions)
                        if chrom_set > 1:
                            # Chromosome Reset
                            max_pos = 0
                            curr_interval = [0, windowsize]
                            yield result_list[:-1]
                            result_list = [result_list, result_list[-1]]
                        elif max_pos >= curr_interval[0] and max_pos < curr_interval[1] and line not in result_list:
                            # If in current interval, do nothing.
                            pass
                        elif max_pos < curr_interval[0]:
                            # If beneath interval - reset to beginning.
                            curr_interval = [0, windowsize]
                        elif max_pos >= curr_interval[1]:
                            # If past interval, iterate and yield result.
                            curr_interval = [curr_interval[1], curr_interval[1] + windowsize]
                            if len(result_list) == 1:
                                yield result_list
                                result_list = []
                            elif len(result_list) > 1:
                                yield result_list[:-1]
                                result_list = [result_list[-1]]
                        line = self.next()

                elif shift_method == ["POS","sliding"]:
                    if len(result_list) == 0:
                        result_list.append(line)
                    while True:
                        positions = [x.POS for x in result_list]
                        max_pos = max(positions)
                        min_pos = min(positions)
                        chrom_num = len(set([o.CHROM for o in result_list])) 
                        if (max_pos - min_pos) > windowsize or chrom_num > 1:
                            result_list.popleft()
                        else:
                            break
                    yield result_list
                    result_list.append(line)

        except StopIteration:
            yield result_list

class variant_set:
    def __init__(self, varset):
        print dir(varset[0])
        print varset[0]
        self.shift_method = shift_method
        self.gt = np.vstack([x.gt_types for x in varset])
        self.segregating_sites = sum([np.any(p) for p in np.equal(3, x)])
        

x = vcf("../test.vcf.gz")
print dir(x)

for i in x.window(windowsize = 100, shift_method = ["SNP","interval"]):
    var = variant_set(i)
    break
    