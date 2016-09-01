from cyvcf2 import VCF as cyvcf2
from cyvcf2 import Variant
from collections import OrderedDict, deque
from itertools import islice
from clint.textui import colored, puts, puts_err, indent
import re
from copy import copy
import os
from . import autoconvert, lev
import numpy as np
import sys
import fileinput
from subprocess import Popen, PIPE, check_output
from reference import resolve_reference_genome
np.set_printoptions(threshold=np.nan)
from collections import defaultdict
from Bio import pairwise2

class cvariant:
    """
    Variant with surrounding region.
    """
    def __init__(self, variant, reference_seq, consensus_seq, template, gt_collection = None):
        self.variant = variant
        self.reference_seq = reference_seq
        self.consensus_seq = consensus_seq
        self.edit_distance = lev(self.reference_seq.lower(), self.consensus_seq.lower())
        if template in [None, "ALT"]:
            self.output = self.consensus_seq
            self.sample = "ALT"
        elif template in ["REF"]:
            self.output = self.reference_seq
            self.sample = "REF"
        else:
            self.output = self.consensus_seq
            self.sample = template
        self.sample_gt = False  # Show sample genotypes.
        if gt_collection:
            self.gt_collection = gt_collection
            self.homozygous_ref = ','.join(self.gt_collection[0])
            self.heterozygous = ','.join(self.gt_collection[1])
            self.homozygous_alt = ','.join(self.gt_collection[3])
            self.missing = ','.join(self.gt_collection[2])

    def __repr__(self):
        CHROM_POS = "{self.variant.CHROM}:{self.variant.POS}".format(**locals())
        return "\t".join([CHROM_POS,
                          self.variant.REF,
                          ','.join(self.variant.ALT),
                          self.sample,
                          str(self.edit_distance),
                          self.output,
                          self.homozygous_ref,
                          self.heterozygous,
                          self.homozygous_alt,
                          self.missing])


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

        self.info_set = [x for x in self.header_iter() if x.type == "INFO"]
        self.filter_set = [x for x in self.header_iter() if x.type == "FILTER"]
        self.format_set = [x for x in self.header_iter() if x.type == "FORMAT"]
        self.header = copy(self.raw_header)


    def insert_header_line(self, header_line):
        header = self.header.splitlines()
        header.insert(len(header)-1, header_line)
        self.header = '\n'.join(header)
        return self.header


    def fetch_variants(self, chrom, start, end):
        """
            Fetch variants by specifying chrom, start, and end.
        """
        region = "{chrom}:{start}-{end}".format(**locals())
        for i in self(region):
            yield i


    def variant_region(self, chrom, start, end, template = None):
        """
            Use samtools to fetch the variant region
            with optional consensus of variants.

            Use template = "ref" to return the reference
            Omit template, or set sample = "alt" to return alternative genotypes.
            Or set template = <sample name> to return a consensus sequence
            for a specific sample.
        """
        sample_flag = ""
        if template == "REF":
            command = "samtools faidx {self.reference_file} {chrom}:{start}-{end}"
        elif template == "ALT" or template is None:
            command = "samtools faidx {self.reference_file} {chrom}:{start}-{end} | bcftools consensus {self.filename}"
        else:
            sample_flag = "--sample=" + template # Get template for sample.
            command = "samtools faidx {self.reference_file} {chrom}:{start}-{end} | bcftools consensus {sample_flag} {self.filename}"
        command = command.format(**locals())
        seq = check_output(command, shell = True)
        return ''.join(seq.splitlines()[1:])


    def fetch_variants_w_consensus(self, chrom, start, end, size = 100, template = None, samples = "ALL"):
        """
            Return a variant with the sequence surrounding it; Both reference
            and alternate or for a set of specified samples.
        """
        if chrom:
            variant_iterator = self.fetch_variants(chrom, start, end)
        else:
            variant_iterator = self

        for variant in variant_iterator:
            start = variant.POS - size
            end = variant.POS + size
            reference_seq = self.variant_region(variant.CHROM,
                                                start,
                                                end,
                                                template="REF")
            if template in ["ALT", "REF"] and samples is None:
                # If template set to REF/ALT and no samples, return template.
                consensus_seq = self.variant_region(variant.CHROM,
                                                    start,
                                                    end,
                                                    template="ALT")
                yield cvariant(variant, reference_seq, consensus_seq, template)
            else:
                # Check that template is sample
                if template not in self.samples + ["REF", "ALT"]:
                    with indent(4):
                        err_msg = "\n" + template + " template not found.\n"
                        puts_err(colored.red(err_msg))
                        exit()
                gt_collection = defaultdict(list)
                if samples != "-":
                    pass

                for vcf_sample, gt in zip(self.samples, variant.gt_types):
                    if vcf_sample in samples:
                        gt_collection[gt].append(vcf_sample)

                consensus_seq = self.variant_region(variant.CHROM,
                                                    variant.POS - size,
                                                    variant.POS + size,
                                                    template=template)
                yield cvariant(variant,
                               reference_seq,
                               consensus_seq,
                               template=template,
                               gt_collection=gt_collection)


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
    """
        Used to modify genotypes by HMM currently...
    """
    def __init__(self, line, sample_names = None):
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
        gt = self.values[index][f_index] = value
        
    def __str__(self):
        self.line = self.line[0:8] + [":".join(self.format_keys)] + [":".join(x) for x in self.values]
        return '\t'.join(self.line)

