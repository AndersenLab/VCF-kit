from cyvcf2 import VCF as cyvcf2
from cyvcf2 import Variant
from collections import OrderedDict, deque
from itertools import islice
from clint.textui import colored, puts, indent
import re
from copy import copy
import os
from . import autoconvert
import numpy as np
import sys
import fileinput
from subprocess import Popen, PIPE, check_output
from reference import resolve_reference_genome
np.set_printoptions(threshold=np.nan)


class vcf(cyvcf2):
    def __init__(self, filename, reference=None):
        if not os.path.isfile(filename) and filename != "-":
            with indent(4):
                exit(puts(colored.red("\nError: " + filename + " does not exist\n")))
        self.filename = filename
        if reference:
            self.reference = reference
            self.reference_file = resolve_reference_genome(reference)

        cyvcf2.__init__(self, self.filename)
        # Check if file exists
        self.n = len(self.samples)  # Number of Samples

        # Meta Data
        comp = re.compile(r'''^##(?P<key>[^<#]+?)=(?P<val>[^<#]+$)''', re.M)
        self.metadata = OrderedDict(comp.findall(self.raw_header))

        # Contigs
        self.contigs = OrderedDict(zip(
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

        self.filter_set = [x for x in self.header_iter() if x.type == "FILTER"]
        self.format_set = [x for x in self.header_iter() if x.type == "FORMAT"]
        self.header = copy(self.raw_header)


    def insert_header_line(self, header_line):
        header = self.header.splitlines()
        header.insert(len(header)-1, header_line)
        self.header = '\n'.join(header)
        return self.header


    def fetch_variants(self, chrom, start, end, samples = ""):
        """
            Use bcftools to fetch variants with view or query 
            from a given interval, optionally subset by sample.
        """
        if samples:
            sample_query = "--samples=" + ','.join(samples)
        region = "{chrom}:{start}-{end}".format(**locals())
        comm = filter(len,["bcftools", "view", "-H", sample_query , self.filename, region])
        comm = Popen(comm, stdout=PIPE, stderr=PIPE)
        for line in comm.stdout:
            yield variant_line(line, samples)


    def consensus(self, chrom, start, end, sample = None):
        """
            Use samtools to fetch consensus; 
        """
        command = "samtools faidx {self.reference_file} {chrom}:{start}-{end} | bcftools consensus {self.filename}"
        command = command.format(**locals())
        seq = check_output(command, shell = True)
        return ''.join(seq.splitlines()[1:])


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
                        if line.POS < result_list.lower_bound:
                            result_list.iterate_interval()
                        elif line.POS >= result_list.upper_bound:
                            yield result_list.filter_within_bounds()
                            result_list.iterate_interval()
                            #for k,i in enumerate(result_list.filter_within_bounds()):
                            #    print k, i.CHROM, i.POS, "\n", result_list.upper_bound
                        if result_list.lower_bound  <= line.POS < result_list.upper_bound:
                            line = self.next()
                            result_list.append(line)
                        last_chrom = line.CHROM
                        #print line.POS
                #
                # POS-Sliding
                #
                elif shift_method == "POS-Sliding":
                    result_list.append(line)
                    result_list.lower_bound = result_list[0].POS - window_size/2
                    if result_list.lower_bound < 0:
                        result_list.lower_bound = 0
                    result_list.upper_bound = result_list.lower_bound + window_size
                    while True:
                        positions = result_list.positions()
                        max_pos = max(positions)
                        min_pos = min(positions)
                        chrom_num = len(set([o.CHROM for o in result_list]))
                        if (max_pos - min_pos) > window_size or chrom_num > 1:
                            result_list.popleft()
                            result_list.iterate_interval() # Updates lower and upper bounds
                        else:
                            break
                    # Update variant interval
                    yield result_list

        except StopIteration:
            if len(result_list) > 0:
                if shift_method == "POS-Interval":
                    yield result_list.filter_within_bounds()
                else:
                    yield result_list

    def output_raw(self):
        """
            Outputs raw vcf
        """
        for hline in self.raw_header.strip("\n").split("\n"):
            yield hline
        for line in self:
            yield(str(line))


class variant_interval(deque):
    def __init__(self, interval = [], window_size = None, step_size = None, shift_method = None, varlist=None, lower_bound = 0):
        self.shift_method = shift_method
        self.window_size = window_size
        self.step_size = step_size
        if self.step_size == None:
            self.step_size = window_size
        self.lower_bound = lower_bound
        self.upper_bound = self.lower_bound + window_size
        if shift_method in ["SNP-Sliding", "SNP-Interval"]:
            super(variant_interval, self).__init__(interval, maxlen = self.window_size)
        else:
            super(variant_interval, self).__init__(interval)

    def positions(self):
        return [x.POS for x in self]

    def iterate_interval(self):
        if self.shift_method == "POS-Sliding":
            mean_pos = np.mean(self.positions())
            half_window = self.window_size/2
            self.lower_bound = int(mean_pos - half_window)
            if self.lower_bound < 0:
                self.lower_bound = 0
            self.upper_bound = int(mean_pos + half_window)
        else:
            self.lower_bound += self.step_size
            self.upper_bound += self.step_size
        self.CHROM = self[0].CHROM
        return self

    def filter_within_bounds(self):
        #remove_list = [x for x in self if x.POS < self.lower_bound]
        #for i in remove_list:
        #    self.remove(i)
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

    #def __repr__(self):
    #    formatted_variants = ["{chrom}:{pos}".format(chrom=self.CHROM, pos=var.POS) for var in self]
    #    formatted_variants = "\t".join(formatted_variants)
    #    interval = self.interval()
    #    return "{0}:{1}-{2}".format(*self.interval()) + " -> " + formatted_variants

class variant_set:
    def __init__(self, ins, samples):
        if type(ins) == str:
            if samples:
                self.variants = [variant_line(x, samples) for x in ins.splitlines()]
            else:
                self.variants =  [variant_line(x) for x in ins.splitlines()]

class variant_line:
    def __init__(self, line, sample_names = None):
        line = str(line).strip().split("\t")
        self.line = line
        self.chrom = line[0]
        self.pos = int(line[1])
        self.ref = line[3]
        self.gt_dict = dict([(n+1,y) for n,y in enumerate(line[4].split(","))] + [(".",".")])
        self.gt_dict[0] = self.ref
        self.alt = self.gt_dict[1]
        self.format_field = line[8].split(":")
        self.has_gt = False
        if "GT" in self.format_field:
            self.gt_loc = self.format_field.index("GT")
            self.has_gt = True
        if sample_names:
            self._sample_to_idx = dict(zip(sample_names, range(9, len(line))))

    def __getitem__(self, i):
        return self.line[i]

    def __setitem__(self, i, val):
        self.line[i] = val

    def append_format_flag(self, flag_name):
        if flag_name not in self.format_field:
            self.format_field += [flag_name]
        return self.format_field.index(flag_name)

    def modify_gt_format(self, i, flag_name, val):
        i += 9
        gt_mod = self.line[i].split(":")
        flag_index = self.append_format_flag(flag_name)
        if len(gt_mod) == flag_index:
            gt_mod += [val]
        else:
            gt_mod[flag_index] = val
        self.line[i] = ':'.join(gt_mod)

    def fetch_gt_from_index(self, i):
        i += 9
        return self.line[i].split(":")[self.gt_loc]

    def gt(self, sample):
        if type(sample) == int:
            return self.fetch_gt_from_index(self._sample_to_idx[sample])
        else:
            return self.line[self._sample_to_idx[sample]].split(":")[self.gt_loc]

    def tgt(self, sample, het = True):
        """
            Return formatted genotype; Homozygous calls only.
        """
        gt = self.line[self._sample_to_idx[sample]].split(":")[self.gt_loc]
        gt_set = re.split("[/|\|]", gt)
        genotype = list(set([self.gt_dict[autoconvert(x)] for x in gt_set]))
        if len(genotype) == 0:
            return genotype[0]
        else:
            return '/'.join(genotype)

    def __getitem__(self, sample):
        return [self.chrom, self.pos, self.tgt(sample)]

    def __repr__(self):
        return self.chrom + ":" + str(self.pos)

    def __str__(self):
        self.line[8] = ':'.join(self.format_field)
        return '\t'.join(self.line)

