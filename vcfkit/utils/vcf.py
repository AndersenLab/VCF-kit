from cyvcf2 import VCF as cyvcf2
from collections import OrderedDict, deque
from itertools import islice
import re
from vcfkit.utils import message
from copy import copy
import os
import numpy as np
from reference import resolve_reference_genome
np.set_printoptions(threshold=np.nan)


class vcf(cyvcf2):
    def __init__(self, filename, reference=None):
        if not os.path.isfile(filename) and filename != "-":
            exit(message("Error: " + filename + " does not exist"))
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
            map(int, re.compile("##contig.*length=([^,>]*?)>").findall(self.raw_header))
        ))

        self.info_set = [x for x in self.header_iter() if x.type == "INFO"]
        self.filter_set = [x for x in self.header_iter() if x.type == "FILTER"]
        self.format_set = [x for x in self.header_iter() if x.type == "FORMAT"]
        self.header = copy(self.raw_header)

    def window(self, shift_method, window_size, step_size=None):
        """
            Generates windows of VCF data by positions or SNP
        """
        result_list = variant_interval(window_size=window_size, step_size=step_size, shift_method=shift_method)
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
                #==============#
                # POS-Interval #
                #==============#
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
                        if result_list.lower_bound <= line.POS < result_list.upper_bound:
                            line = self.next()
                            result_list.append(line)
                        last_chrom = line.CHROM
                #=============#
                # POS-Sliding #
                #=============#
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
                            result_list.iterate_interval()  # Updates lower and upper bounds
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
    def __init__(self, interval=[], window_size=None, step_size=None, shift_method=None, varlist=None, lower_bound=0):
        self.shift_method = shift_method
        self.window_size = window_size
        self.step_size = step_size
        if self.step_size is None:
            self.step_size = window_size
        self.lower_bound = lower_bound
        self.upper_bound = self.lower_bound + window_size
        if shift_method in ["SNP-Sliding", "SNP-Interval"]:
            super(variant_interval, self).__init__(interval, maxlen=self.window_size)
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
        cut = [x for x in self if x.POS >= self.lower_bound and x.POS < self.upper_bound]
        return variant_interval(interval=cut,
                                window_size=self.window_size,
                                shift_method=self.shift_method,
                                lower_bound=self.lower_bound,
                                step_size=self.step_size)

    def __getitem__(self, index):
        if isinstance(index, slice):
            cut = islice(self, index.start, index.stop, index.step)
            return variant_interval(interval=cut,
                                    window_size=self.window_size,
                                    shift_method=self.shift_method,
                                    lower_bound=self.lower_bound,
                                    step_size=self.step_size)
        return deque.__getitem__(self, index)

    def get_last(self):
        """
            Removes all items from the variant interval and
            sets to last item
        """
        super(variant_interval, self).__init__(self[len(self)-1:len(self)], maxlen=self.window_size)

    def unique_chroms(self):
        """
            Checks to see whether the variants
            all belong to the same chromosome
        """
        return len(set([var.CHROM for var in self])) == 1 and len(self) > 0


class variant_set:
    def __init__(self, ins, samples):
        if type(ins) == str:
            if samples:
                self.variants = [variant_line(x, samples) for x in ins.splitlines()]
            else:
                self.variants = [variant_line(x) for x in ins.splitlines()]


class variant_line:
    """
        Used to modify genotypes by HMM currently...
    """
    def __init__(self, line, sample_names=None):
        self.line = str(line).strip().split("\t")
        self.chrom = self.line[0]
        self.pos = int(self.line[1])
        self.chrom_pos = (self.chrom, self.pos)
        self.format_keys = self.line[8].split(":")
        self.values = [x.split(":") for x in self.line[9:]]

    def get_gt(self, format_field, index):
        f_index = self.format_keys.index(format_field)
        return self.values[index][f_index]

    def set_gt(self, format_field, index, value):
        if format_field not in self.format_keys:
            self.format_keys += [format_field]
            self.values = [x + ["."] for x in self.values]
        f_index = self.format_keys.index(format_field)
        self.values[index][f_index] = value

    def __str__(self):
        self.line = self.line[0:8] + [":".join(self.format_keys)] + [":".join(x) for x in self.values]
        return '\t'.join(self.line)
