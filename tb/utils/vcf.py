from cyvcf2 import VCF as cyvcf2
from collections import OrderedDict, deque
from itertools import islice
from clint.textui import colored, puts, indent
import re
from copy import copy
import os
import numpy as np
import sys
import fileinput
from subprocess import Popen, PIPE
np.set_printoptions(threshold=np.nan)


class vcf(cyvcf2):
    def __init__(self, filename, reference=None):
        if not os.path.isfile(filename) and filename != "-":
            with indent(4):
                exit(puts(colored.red("\nError: " + filename + " does not exist\n")))
        self.filename = filename
        if reference:
            self.reference = reference

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

        self.header = copy(self.raw_header)

    def insert_header_line(self, header_line):
        header = self.header.splitlines()
        header.insert(len(header)-1, header_line)
        self.header = '\n'.join(header)
        return self.header


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

    def fasta(self, region, sample = "", het=False, hom=True):
        """
            Fetch a fasta sequence; spike with variant calls for a given sample
        """
        chrom, start, end = re.split("[:-]", region)
        start, end = int(start), int(end)
        fasta, err = Popen(["samtools", "faidx", self.reference, region], stdout = PIPE, stderr = PIPE).communicate()
        fasta = ''.join(fasta.split("\n")[1:])
        var_col = 3
        comm = ["bcftools", "view", "-H", self.filename, region]

        # Use het?
        use_het = 3
        if het:
            use_het = 4

        ref_alt = {"./.": 3, "0/0": 3, "0/1": use_het, "1/1": 4}
        if sample:
            if sample != "ALT":
                comm.insert(2,"--samples=" + sample)
            var_col = 4
        if sample:
            variants, err = Popen(comm, stdout = PIPE, stderr = PIPE).communicate()
            variants = [x.split("\t") for x in variants.strip().split("\n")]
            if variants[0][0]:
                for variant in variants:
                    chrom, pos = variant[0], int(variant[1])
                    fm = variant[8]
                    gt = variant[9]
                    GT_loc = fm.split(":").index("GT")
                    gt = gt.split(":")[GT_loc]
                    if gt in ref_alt and sample != "ALT":
                        var_col = ref_alt[gt]
                    pos = pos - start
                    fasta = fasta[:pos] + variant[var_col] + fasta[pos+1:]
                return fasta
        else:
            # Return reference sequence
            return fasta


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

class variant_line:
    def __init__(self, line):
        line = str(line).strip().split("\t")
        self.line = line
        self.chrom = line[0]
        self.pos = int(line[1])
        self.format_field = line[8].split(":")
        self.has_gt = False
        if "GT" in self.format_field:
            self.gt_loc = self.format_field.index("GT")
            self.has_gt = True

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

    def __str__(self):
        self.line[8] = ':'.join(self.format_field)
        return '\t'.join(self.line)

