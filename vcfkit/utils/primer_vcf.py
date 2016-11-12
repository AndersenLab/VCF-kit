from cyvcf2 import VCF as cyvcf2
from clint.textui import colored, puts, indent
import re
from collections import OrderedDict, defaultdict
from vcfkit.utils import lev, message
from copy import copy
import os
import numpy as np
from utils.primer3 import primer3
from subprocess import Popen, PIPE, check_output
from reference import resolve_reference_genome
np.set_printoptions(threshold=np.nan)

debug = None


class cvariant:
    """
    Variant with surrounding region.
    ref_seq = Reference Sequence of variant
    alt_seq = Alternative Sequence of variant (as determined by template)
    box_seq = reference sequence with alternative variant boxed.
    """
    def __init__(self, variant, ref_seq, alt_seq, box_seq, template, gt_collection=None):
        self.variant = variant
        self.ref_seq = ref_seq
        self.alt_seq = alt_seq
        self.box_seq = box_seq
        self.edit_distance = lev(self.ref_seq.lower(), self.alt_seq.lower())
        if template in [None, "ALT"]:
            self.output = self.alt_seq
            self.sample = "ALT"
        elif template in ["REF"]:
            self.output = self.ref_seq
            self.sample = "REF"
        else:
            self.output = self.alt_seq
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
        output = [CHROM_POS,
                  self.variant.REF,
                  ','.join(self.variant.ALT),
                  self.sample,
                  str(self.edit_distance),
                  self.output,
                  self.homozygous_ref,
                  self.heterozygous,
                  self.homozygous_alt,
                  self.missing]
        if self.box_seq:
            output.insert(6, self.box_seq)
        return "\t".join(output)


class primer_vcf(cyvcf2):
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
            map(int, re.compile("##contig.*length=(.*?)>").findall(self.raw_header))
        ))

        self.info_set = [x for x in self.header_iter() if x.type == "INFO"]
        self.filter_set = [x for x in self.header_iter() if x.type == "FILTER"]
        self.format_set = [x for x in self.header_iter() if x.type == "FORMAT"]
        self.header = copy(self.raw_header)

    def variant_iterator(self):
        """
            Iterate over variants in region is specified;
            Otherwise start at the beginning.
        """
        if self.region:
            var_region = self.fetch_variants(self.region)
        else:
            var_region = self

        for variant in var_region:
            size = self.size
            start = variant.POS - size
            end = variant.POS + size

            if self.template not in self.samples + ["REF", "ALT"]:
                exit(message(self.template + " template not found."))
            gt_collection = defaultdict(list)

            for vcf_sample, gt in zip(self.samples, variant.gt_types):
                if vcf_sample in self.output_samples:
                    gt_collection[gt].append(vcf_sample)

            # Retrieve ref and alt sequences; Define box_seq as none.
            ref_seq, alt_seq = self.variant_region(variant, size)
            box_seq = None

            if self.box_variants:
                a, b = alt_seq[:size], alt_seq[size+len(variant.REF):]
                alt_variants = '/'.join(variant.ALT)
                box_seq = "{a}[{variant.REF}/{alt_variants}]{b}".format(**locals())

            yield cvariant(variant,
                           ref_seq,
                           alt_seq,
                           box_seq,
                           template=self.template,
                           gt_collection=gt_collection)

    def fetch_variants(self, region=None):
        """
            Fetch variants using the specified region
        """
        if region is None:
            region = self.region
        for i in self(region):
            yield i

    def variant_region(self, variant, size):
        """
            Use samtools to fetch the variant region
            with optional consensus of variants.

            Use template = "ref" to return the reference
            Omit template, or set sample = "alt" to return alternative genotypes.
            Or set template = <sample name> to return a consensus sequence
            for a specific sample.
        """
        chrom = variant.CHROM
        start = variant.POS - size
        end = variant.POS + size
        sample_flag = ""
        template = self.template

        # Make sure index exists
        if not os.path.isfile(self.filename + ".csi"):
            message("VCF is not indexed; Indexing.")
            comm = "bcftools index {f}".format(f=self.filename)
            out = check_output(comm, shell=True)
            message(out)

        # Retreive the reference sequence
        ref_seq_command = "samtools faidx {self.reference_file} {chrom}:{start}-{end}"
        ref_seq_command = ref_seq_command.format(**locals())
        ref_seq = check_output(ref_seq_command, shell=True)
        ref_seq = ''.join(ref_seq.splitlines()[1:])

        if template == "REF":
            command = "samtools faidx {self.reference_file} {chrom}:{start}-{end}"
        elif template == "ALT" or template is None:
            command = "samtools faidx {self.reference_file} {chrom}:{start}-{end} | bcftools consensus {self.filename}"
        else:
            sample_flag = "--sample=" + template  # Get template for sample.
            command = "samtools faidx {self.reference_file} {chrom}:{start}-{end} | bcftools consensus {sample_flag} {self.filename}"
        command = command.format(**locals())
        alt_seq = check_output(command, shell=True)
        alt_seq = ''.join(alt_seq.splitlines()[1:])
        return ref_seq, alt_seq

    def print_primer_header(self):
        header = ["CHROM_POS",
                  "REF",
                  "ALT",
                  "SAMPLE",
                  "EDIT_DISTANCE",
                  "SEQUENCE",
                  "0/0",
                  "0/1",
                  "1/1",
                  "./."]
        if self.box_variants:
            header.insert(6, "SEQUENCE_BOX")
        if self.mode == "template":
            print('\t'.join(header))


    def fetch_template(self):
        """
            Return a variant with the sequence surrounding it; Both reference
            and alternate or for a set of specified samples.

            sample_output: List of samples to output w/ their genotypes.
        """
        for variant in self.variant_iterator():
            yield variant


    def fetch_indel_primers(self):
        for variant in self.variant_iterator():
            p3 = primer3()
            yield p3.fetch_primers(variant.ref_seq)
            # Take processed variant and feed to indel processor.

