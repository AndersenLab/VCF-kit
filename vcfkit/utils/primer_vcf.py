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
import signal
signal.signal(signal.SIGINT, lambda x,y: sys.exit(0))


high_fidelity = ["AgeI",
                 "ApoI",
                 "BamHI",
                 "BbsI",
                 "BmtI",
                 "BsaI",
                 "BsiWI",
                 "BsrGI",
                 "BstEII",
                 "BstZ17I",
                 "DraIII",
                 "EagI",
                 "EcoRI",
                 "EcoRV",
                 "HindIII",
                 "KpnI",
                 "MfeI",
                 "MluI",
                 "NcoI",
                 "NheI",
                 "NotI",
                 "NruI",
                 "NsiI",
                 "PstI",
                 "PvuI",
                 "PvuII",
                 "SacI",
                 "SalI",
                 "SbfI",
                 "ScaI",
                 "SpeI",
                 "SphI",
                 "SspI",
                 "StyI"]

debug = None


class cvariant:
    """
    Variant with surrounding region.
    ref_seq = Reference Sequence of variant
    alt_seq = Alternative Sequence of variant (as determined by template)
    box_seq = reference sequence with alternative variant boxed.
    """
    def __init__(self, variant):
        for i in filter(lambda x: x.startswith("_") is False, dir(variant)):
            setattr(self, i, getattr(variant, i))
      
    def __repr__(self):
        CHROM_POS = "{self.CHROM}:{self.POS}".format(**locals())

        if self.template in [None, "ALT"]:
            self.output = self.alt_seq
            self.sample = "ALT"
        elif self.template in ["REF"]:
            self.output = self.ref_seq
            self.sample = "REF"
        else:
            self.output = self.alt_seq
            self.sample = self.template

        output = [CHROM_POS,
                  self.REF,
                  ','.join(self.ALT),
                  self.sample,
                  str(self.edit_distance),
                  self.output.tostring(),
                  self.homozygous_ref,
                  self.heterozygous,
                  self.homozygous_alt,
                  self.missing]
        if self.box_seq:
            output.insert(6, self.box_seq)
        return "\t".join(output)



from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as DNA_SET
from Bio.Restriction import AllEnzymes

class restriction_sites:

    """
        Container for determining restriction sites
        and modifying their coordinates
    """

    def __init__(self, ref_seq, alt_seq):
        self.ref_sites = dict(AllEnzymes.search(self.ref_seq).items())
        self.alt_sites = dict(AllEnzymes.search(self.alt_seq).items())
        self.ref_sites_diff = {k: v for k, v in self.ref_sites.items() if len(v) > 0 and len(v) <= 3 and
                               abs(len(self.ref_sites[k]) - len(self.alt_sites[k])) == 1}
        self.alt_sites_diff = {k: self.alt_sites[k] for k, v in self.ref_sites_diff.items()}

    def shift_positions(self, shift):
        for k in self.ref_sites_diff.keys():
            for n, s1 in enumerate(self.ref_sites_diff[k]):
                self.ref_sites_diff[k][n] = s1 + shift
            for n, s2 in enumerate(self.alt_sites_diff[k]):
                self.alt_sites_diff[k][n] = s2 + shift

    def __len__(self):
        return len(self.ref_sites_diff.keys())


    def extract_restriction(self, window=2000):
        for varset in self.window(window_size=window, step_size=1, shift_method="POS-Sliding"):
            start = varset.lower_bound
            end = varset.upper_bound
            CHROM = varset[0].CHROM
            ref_seq = self.reference[CHROM][start:end].seq
            alt_seq = ref_seq
            for var in varset:
                # Spike all alt variants within sequence.
                alt_seq = alt_seq[:var.POS - start] + var.ALT[0] + alt_seq[var.POS - start + 1:]
            rsites = restriction_sites(ref_seq, alt_seq)
            if len(rsites) > 0:
                primer_set = primer3()
                for primer in primer_set.fetch_primers(ref_seq):
                    yield primer, restriction_sites, start, end




class primer_vcf(cyvcf2):
    def __init__(self, filename, reference=None, template="ALT"):
        if not os.path.isfile(filename) and filename != "-":
            exit(message("Error: " + filename + " does not exist"))
        self.filename = filename
        self.template = template
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

        # Setup template
        if self.template not in self.samples + ["REF", "ALT"]:
            exit(message(self.template + " template not found."))


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
            nz = cvariant(variant)
            nz.template = self.template
            nz.mode = self.mode

            gt_collection = defaultdict(list)

            for vcf_sample, gt in zip(self.samples, variant.gt_types):
                if vcf_sample in self.output_samples:
                    gt_collection[gt].append(vcf_sample)

            # Determine region size
            nz.indel_size = 0
            nz.region_size = self.region_size
            if nz.is_indel:
                nz.indel_size =  abs(len(nz.REF) - len(nz.ALT[0]))
                nz.region_size = nz.indel_size*2 + 150


            # Retrieve ref and alt sequences; Define box_seq as none.
            ref_seq, alt_seq = self.variant_region(nz)
            nz.ref_seq = ref_seq
            nz.alt_seq = alt_seq
            nz.edit_distance = lev(ref_seq.tostring().lower(), alt_seq.tostring().lower())

            if gt_collection:
                nz.gt_collection = gt_collection
                nz.homozygous_ref = ','.join(nz.gt_collection[0])
                nz.heterozygous = ','.join(nz.gt_collection[1])
                nz.homozygous_alt = ','.join(nz.gt_collection[3])
                nz.missing = ','.join(nz.gt_collection[2])

            nz.box_seq = None
            if self.box_variants:
                a, b = alt_seq[:nz.region_size], alt_seq[nz.region_size+len(variant.REF):]
                alt_variants = '/'.join(variant.ALT)
                nz.box_seq = "{a}[{variant.REF}/{alt_variants}]{b}".format(**locals())

            yield nz

    def fetch_variants(self, region=None):
        """
            Fetch variants using the specified region
        """
        if region is None:
            region = self.region
        for i in self(region):
            yield i

    def variant_region(self, variant):
        """
            Use samtools to fetch the variant region
            with optional consensus of variants.

            Use template = "ref" to return the reference
            Omit template, or set sample = "alt" to return alternative genotypes.
            Or set template = <sample name> to return a consensus sequence
            for a specific sample.
        """
        chrom = variant.CHROM
        start = variant.POS - variant.region_size
        end = variant.POS + variant.region_size
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
        ref_seq = Seq(''.join(ref_seq.splitlines()[1:]), DNA_SET)

        if template == "REF":
            command = "samtools faidx {self.reference_file} {chrom}:{start}-{end}"
        elif template == "ALT" or template is None:
            command = "samtools faidx {self.reference_file} {chrom}:{start}-{end} | bcftools consensus {self.filename}"
        else:
            sample_flag = "--sample=" + template  # Get template for sample.
            command = "samtools faidx {self.reference_file} {chrom}:{start}-{end} | bcftools consensus {sample_flag} {self.filename}"
        command = command.format(**locals())
        try:
            alt_seq = check_output(command, shell=True)
            alt_seq = Seq(''.join(alt_seq.splitlines()[1:]), DNA_SET)
        except:
            alt_seq = ''
        return ref_seq, alt_seq

    def print_primer_header(self):
        header = ["CHROM_POS",
                  "REF",
                  "ALT",
                  "SAMPLE",
                  "EDIT_DISTANCE",
                  "TEMPLATE",
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
        p3 = primer3()
        for variant in self.variant_iterator():
            if 30 <= variant.indel_size <= 500:
                low = variant.region_size/2
                high = variant.region_size
                p3.PRIMER_PRODUCT_SIZE_RANGE = "{low}-{high}".format(**locals())
                print list(p3.fetch_primers(variant.ref_seq))
                yield p3.fetch_primers(variant.ref_seq)
                print(variant)
                # Take processed variant and feed to indel processor.

    def fetch_snpsnp_primers(self):
        for variant in self.variant_iterator():
            yield variant

    def fetch_sanger_primers(self):
        for variant in self.variant_iterator():
            yield variant
