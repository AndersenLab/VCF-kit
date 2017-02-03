from cyvcf2 import VCF as cyvcf2
from clint.textui import colored, puts, indent
import re
from collections import OrderedDict, defaultdict
from vcfkit.utils import lev, message
from copy import copy
import os
import numpy as np
from cyvcf2 import VCFReader
from vcfkit.utils.primer3 import primer3
from subprocess import Popen, PIPE, check_output
from reference import resolve_reference_genome
np.set_printoptions(threshold=np.nan)
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as DNA_SET
from Bio.Restriction import AllEnzymes, CommOnly, RestrictionBatch

# Global flag for header output
header_printed = False

high_fidelity = RestrictionBatch(["AgeI",
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
                 "StyI"])

debug = None


class cvariant:
    """
        Mutable variant object
    """
    def __init__(self, variant):
        for i in filter(lambda x: x.startswith("_") is False, dir(variant)):
            setattr(self, i, getattr(variant, i))


class template:

    def __init__(self, variant, reference_file, vcf, use_template, nprimers):
        # Copy variant attributes over
        for i in filter(lambda x: x.startswith("_") is False, dir(variant)):
            setattr(self, i, getattr(variant, i))
        self.output_samples = vcf.output_samples
        self.reference_file = reference_file
        self.region_start = variant.POS - variant.region_size
        if self.region_start <= 0:
            self.region_start = 1
        self.region_end = variant.POS + variant.region_size - 1
        self.region = "{self.CHROM}:{self.region_start}-{self.region_end}".format(**locals())
        
        # Retrieve sequences
        self.ref_seq = self.fetch_sequence("REF")
        if use_template == 'REF':
            self.alt_seq = self.ref_seq
        else:
            self.alt_seq = self.fetch_sequence(use_template)
        self.edit_distance = lev(str(self.ref_seq).lower(), str(self.alt_seq).lower())

        # Retrieve primers
        if self.mode != 'template':
            p3 = primer3(reference_file)
            p3.PRIMER_NUM_RETURN = nprimers
            # Modify amplicon size if indel.
            if self.mode == 'indel':
                max_amplicon_length = max([len(self.ref_seq), len(self.alt_seq)])
                min_amplicon_length = max_amplicon_length - 200
                if (min_amplicon_length < 50):
                    min_amplicon_length = 50
                p3.PRIMER_PRODUCT_SIZE_RANGE = "{min_amplicon_length}-{max_amplicon_length}".format(**locals())
            elif self.mode == 'sanger':
                p3.PRIMER_PRODUCT_SIZE_RANGE = "{self.amplicon_lower}-{self.amplicon_upper}".format(**locals())
            primers = p3.fetch_primers(self.ref_seq, self.CHROM, self.region_start)
            if primers:
                self.primers = [x for x in primers if x.filter_primer_group()]
            else:
                self.primers = []

        # Identify restriction site locations
        if self.mode == 'snip':
            ref_seq = str(self.ref_seq)
            # Generate a primary variant only seq to identify restriction
            # sites that target ALT sites.
            self.primary_variant_seq = Seq(ref_seq[0:500] + self.ALT[0] + ref_seq[501:])
            self.fetch_restriction_sites(vcf.enzymes)


    def fetch_sequence(self, use_template):
        """
            Fetches sequence surrounding variant for REF, ALT, or given sample.
        """
        sample_flag = ""
        if use_template == "REF":
            command = "samtools faidx {self.reference_file} {self.region}"
        elif use_template == "ALT" or use_template is None:
            command = "samtools faidx {self.reference_file} {self.region} | bcftools consensus {self.filename}"
        else:
            sample_flag = "--sample=" + use_template  # Get use_template for sample.
            command = "samtools faidx {self.reference_file} {self.region} | bcftools consensus {sample_flag} {self.filename}"
        command = command.format(**locals())
        try:
            seq = check_output(command, shell=True)
            seq = Seq(''.join(seq.splitlines()[1:]), DNA_SET)
        except:
            seq = Seq('')
        return seq


    def fetch_restriction_sites(self, enzymes = "Common"):
        """ 
            Spike in target variant first, generate list 
            restriction enzymes that will work.
        """
        if enzymes == "ALL":
            enzyme_group = AllEnzymes
        elif enzymes == "Common":
            enzyme_group = CommOnly
        elif enzymes == "HF":
            enzyme_group = high_fidelity
        else:
            enzyme_group = RestrictionBatch(enzymes.split(","))

        # Filter ambiguous cutters
        enzyme_group = RestrictionBatch([x for x in enzyme_group if x.is_ambiguous() is False])

        # Calculate rflps for ALT sites only
        self.ref_sites = dict(enzyme_group.search(self.ref_seq).items())
        self.primary_variant_sites = dict(enzyme_group.search(self.primary_variant_seq).items())
        self.rflps = {k: (self.ref_sites[k], self.primary_variant_sites[k]) for k, v in self.ref_sites.items() 
                        if len(v) > 0 and len(v) <= 3
                        and self.ref_sites[k] != self.primary_variant_sites[k]}


    def calculate_cuts(self, length, positions):
        """
            Given a template and cut sites,
            returns resulting fragments and their size.
        """
        if positions:
            init, cuts = 0, []
            for cut in positions:
                cut_len = cut - init
                cuts.append(cut_len)
                init = init + cut_len
            # Add last piece in
            cuts.append(length - positions[-1])
            return cuts
        else:
            return [length]


    def fetch_variant_count(self, region, sample_set = None):
        """
            Fetches variant count in a region.
        """
        v = VCFReader(self.filename)
        sample_indices = [v.samples.index(x) for x in sample_set]
        vc = 0
        for i in v(self.region):
            gt_len = len(set([i.gt_types[x] for x in sample_indices]))
            if gt_len > 1:
                vc += 1
        return vc



    def print_header(self, hout):
        global header_printed
        if header_printed is False:
            print('\t'.join(map(str, hout)))
            header_printed = True

    def out(self):

        hout = ["CHROM",
                "POS",
                "region",
                "REF",
                "ALT",
                "template_sample"]

        out = [self.CHROM,
               self.POS,
               self.region,
               self.REF,
               ','.join(self.ALT),
               self.use_template]
        homozygous_ref = ','.join(map(str,self.gt_collection[0]))
        homozygous_alt = ','.join(map(str,self.gt_collection[3]))

        if self.mode == 'snip':
            hout += ["variant_count",
                     "ref_sites",
                     "alt_sites",
                     "restriction_enzyme",
                     "restriction_site",
                     "restriction_site_length",
                     "primer_left",
                     "primer_right",
                     "melting_temperature",
                     "amplicon_length",
                     "amplicon_region",
                     "amplicon_sequence",
                     "0/0",
                     "1/1",
                     "polymorphic"]

            self.print_header(hout)
            for primer_group in self.primers:
                amp_start = self.region_start + primer_group.primer_left.START
                amp_end = self.region_start + primer_group.primer_right.END
                amplicon_length = len(primer_group.amplicon)
                primer_group.amplicon_region = "{self.CHROM}:{amp_start}-{amp_end}".format(**locals())
                if primer_group is not None:
                    for rflp, cut_sites in self.rflps.items():
                        # Shift cut sites
                        ref_cut_sites = [x - primer_group.primer_left.START for x in cut_sites[0]]
                        alt_cut_sites = [x - primer_group.primer_left.START for x in cut_sites[1]]
                        # Calculate Cuts and shift positions
                        ref_product_sizes = self.calculate_cuts(amplicon_length, ref_cut_sites)
                        alt_product_sizes = self.calculate_cuts(amplicon_length, alt_cut_sites)

                        # FILTER: If product sizes are less than 25, ignore.
                        product_sizes = ref_product_sizes + alt_product_sizes
                        if any([x < 25 for x in product_sizes]):
                            break

                        # FILTER: Remove sites where the number of restriction sites does not change.
                        if len(ref_cut_sites) == len(alt_cut_sites):
                            break

                        ref_cut_sites = ','.join(map(str, ref_cut_sites))
                        ref_product_sizes = ','.join(map(str, ref_product_sizes))
                        alt_cut_sites = ','.join(map(str, alt_cut_sites))
                        alt_product_sizes = ','.join(map(str, alt_product_sizes))

                        # Fetch variant Count
                        variant_count = self.fetch_variant_count(primer_group.amplicon_region, self.output_samples)

                        # out_primers represents variables that are specific to the primer/rflp
                        out_primers = out + [variant_count,
                                             ref_cut_sites + ":" + ref_product_sizes,
                                             alt_cut_sites + ":" + alt_product_sizes,
                                             rflp,
                                             rflp.site,
                                             rflp.size,
                                             primer_group.primer_left,
                                             primer_group.primer_right,
                                             primer_group.primer_tm,
                                             primer_group.amplicon_length,
                                             primer_group.amplicon_region,
                                             primer_group.amplicon,
                                             homozygous_ref,
                                             homozygous_alt,
                                             self.is_polymorphic]
                        print('\t'.join(map(str, out_primers)))
        elif self.mode == 'indel':
            hout += ["variant_count",
                     "indel_size",
                     "indel_type",
                     "primer_left",
                     "primer_right",
                     "melting_temperature",
                     "REF_product_size",
                     "ALT_product_size",
                     "amplicon_region",
                     "amplicon_sequence",
                     "0/0",
                     "1/1",
                     "polymorphic"]

            self.print_header(hout)
            for primer_group in self.primers:
                # Fetch variant Count
                variant_count = self.fetch_variant_count(primer_group.amplicon_region, self.output_samples)
                out_indel = out + [variant_count,
                                   self.indel_size,
                                   self.indel_type,
                                   primer_group.primer_left,
                                   primer_group.primer_right,
                                   primer_group.primer_tm,
                                   primer_group.amplicon_length,
                                   primer_group.amplicon_length - self.indel_size,
                                   primer_group.amplicon_region,
                                   primer_group.amplicon,
                                   homozygous_ref,
                                   homozygous_alt,
                                   self.is_polymorphic]
                print('\t'.join(map(str, out_indel)))

        elif self.mode == 'sanger':
            hout += ["variant_count",
                     "variant_distance",
                     "primer_left",
                     "primer_right",
                     "melting_temperature",
                     "amplicon_length",
                     "amplicon_region",
                     "amplicon_sequence",
                     "0/0",
                     "1/1",
                     "polymorphic"]
            self.print_header(hout)
            for primer_group in self.primers:
                amp_start = self.region_start + primer_group.primer_left.START
                # Fetch variant Count
                variant_count = self.fetch_variant_count(primer_group.amplicon_region, self.output_samples)
                out_indel = out + [variant_count,
                                   self.POS - amp_start,
                                   primer_group.primer_left,
                                   primer_group.primer_right,
                                   primer_group.primer_tm,
                                   primer_group.amplicon_length,
                                   primer_group.amplicon_region,
                                   primer_group.amplicon,
                                   homozygous_ref,
                                   homozygous_alt,
                                   self.is_polymorphic]
                print('\t'.join(map(str, out_indel)))

        else:
            # Output template
            hout += ["template"]
            out += [self.alt_seq]
            self.print_header(hout)
            print('\t'.join(map(str, out)))


class primer_vcf(cyvcf2):
    def __init__(self, filename, reference=None, use_template="ALT", polymorphic=True):
        if not os.path.isfile(filename) and filename != "-":
            exit(message("Error: " + filename + " does not exist"))
        self.filename = filename
        self.use_template = use_template
        self.use_polymorphic = polymorphic
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

        # Make sure index exists for this VCF.
        if not os.path.isfile(self.filename + ".csi"):
            message("VCF is not indexed; Indexing.")
            comm = "bcftools index {f}".format(f=self.filename)
            out = check_output(comm, shell=True)
            message(out)

        # Setup use_template
        if self.use_template not in self.samples + ["REF", "ALT"]:
            exit(message(self.use_template + " use_template not found."))


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

            # Skip variants that are not applicable
            if self.mode == 'indel':
                if variant.is_indel is False:
                    continue
            elif self.mode == 'snip':
                if variant.is_snp is False:
                    continue

            nz = cvariant(variant)

            # Process Indels
            if nz.is_indel:
                # Size is dynamically set for indels
                nz.indel_size = len(nz.REF) - len(nz.ALT[0])
                nz.indel_size_abs = abs(len(nz.ALT[0]) - len(nz.REF))
                # Indel type:
                if nz.indel_size > 0:
                    nz.indel_type = "deletion"
                else:
                    nz.indel_type = "insertion"

            if self.mode == 'snip':
                nz.region_size = 500
            elif self.mode == 'indel':
                nz.region_size = nz.indel_size * 2 + 150
                if nz.indel_size_abs <= 25:
                    continue
            elif self.mode in ['sanger', 'template']:
                nz.region_size = self.region_size
            nz.use_template = self.use_template
            nz.mode = self.mode
            nz.filename = self.filename
            nz.use_polymorphic = self.use_polymorphic

            gt_collection = defaultdict(list)

            for vcf_sample, gt in zip(self.samples, variant.gt_types):
                if vcf_sample in self.output_samples:
                    gt_collection[gt].append(vcf_sample)

            if gt_collection:
                nz.gt_collection = gt_collection
                nz.homozygous_ref = nz.gt_collection[0]
                nz.heterozygous = nz.gt_collection[1]
                nz.homozygous_alt = nz.gt_collection[3]
                nz.missing = nz.gt_collection[2]
            is_polymorphic = all(map(lambda x: len(x) > 0,[gt_collection[0], gt_collection[3]]))

            nz.is_polymorphic = is_polymorphic
            if self.use_polymorphic and not is_polymorphic:
                continue

            if self.mode == 'sanger':
                nz.amplicon_lower = self.amplicon_lower
                nz.amplicon_upper = self.amplicon_upper

            t = template(nz, self.reference_file, self, self.use_template, self.nprimers)
            yield t

    def fetch_variants(self, region=None):
        """
            Fetch variants using the specified region
        """
        if region is None:
            region = self.region
        for i in self(region):
            yield i

