#! /usr/bin/env python
"""
usage:
  vk primer template [options] <vcf> 
  vk primer sanger   [options] <vcf>
  vk primer snipsnp  [options] <vcf>
  vk primer indel    [options] <vcf>

Example

options:
  -h --help                   Show this screen.
  --version                   Show version.
  --ref=<reference>           Reference Genome
  --region=<region>           Restrict to region.
  --samples=<samples>         Output genotypes for a sample or set of samples. [default: ALL]
  --template=<sample>         Sequence to use for template. Use REF, ALT, or sample name. [default: ALT]
  --mark-variant              Use brackets to indicate the location of the variant in output.
  --size=<int>                Amplicon size/2 (Upstream and downstream length) [default: 50]


"""
from docopt import docopt
from clint.textui import colored, puts, indent
from utils import parse_region, message
from utils.vcf import *
from utils.reference import *
from utils.fasta import *
import sys
from utils.primer3 import primer3
import os
from glob import glob
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as DNA_SET
from Bio.Restriction import AllEnzymes
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)


debug = None
if len(sys.argv) == 1:
    debug = ['primer', "--ref=WBcel235", "test.vcf.gz"]


class restriction_sites:
    """
        Container for determining restriction sites
        and modifying their coordinates
    """
    def __init__(self, ref_seq, alt_seq):
        self.ref_seq = Seq(ref_seq, DNA_SET)
        self.alt_seq = Seq(alt_seq, DNA_SET)
        self.ref_sites = dict(AllEnzymes.search(self.ref_seq).items())
        self.alt_sites = dict(AllEnzymes.search(self.alt_seq).items())
        self.ref_sites_diff = {k:v for k,v in self.ref_sites.items() if len(v) > 0 and len(v) <= 3 and
                                                         abs(len(self.ref_sites[k]) - len(self.alt_sites[k])) == 1}
        self.alt_sites_diff = {k:self.alt_sites[k] for k,v in self.ref_sites_diff.items()}

    def shift_positions(self, shift):
        for k in self.ref_sites_diff.keys():
            for n, s1 in enumerate(self.ref_sites_diff[k]):
                self.ref_sites_diff[k][n] = s1 + shift
            for n, s2 in enumerate(self.alt_sites_diff[k]):
                self.alt_sites_diff[k][n] = s2 + shift

    def __len__(self):
        return len(self.ref_sites_diff.keys())


# class seq_vcf(vcf):
#     """
#         Subclass of vcf that extracts genomic context of variants from reference.
#     """
#     def __init__(self, filename, reference):
#         vcf.__init__(self,filename)
#         if not os.path.isfile(reference):
#             module_path = os.path.split(os.path.realpath(__file__))[0]
#             self.reference = module_path + "/reference/" + args["--ref"] + ".fa.gz"
#             if not os.path.isfile(self.reference):
#                 with indent(4):
#                     puts(colored.red("\nGenome does not exist!\n"))
#                     puts(colored.blue("Available Genomes:\n"))
#                 with indent(8):
#                     for i in glob(module_path + "/reference/*.fa.gz"):
#                         genome = os.path.split(i)[1].replace(".fa.gz", "")
#                         puts(genome)
#                 with indent(4):
#                     puts(colored.blue("\nDownload more with vk genome\n"))
#                     exit()
#             self.reference = Fasta(self.reference)
#             if len([x for x in self.contigs.keys() if x not in self.reference.keys()]) > 0:
#                 with indent(4):
#                     puts(colored.yellow("""\nWarning: VCF / Reference contigs don't match\nUsing contig names within VCF.\n"""))
#                     self.reference.alt_contig_names = dict(zip(self.contigs.keys(), self.reference.keys(),))


    def extract_restriction(self, window = 2000):
        for varset in self.window(window_size=window, step_size = 1, shift_method="POS-Sliding"):
            start = varset.lower_bound
            end = varset.upper_bound
            CHROM = varset[0].CHROM
            ref_seq = self.reference[CHROM][start:end].seq
            alt_seq = ref_seq
            for var in varset:
                # Spike all alt variants within sequence.
                alt_seq = alt_seq[:var.POS-start] + var.ALT[0] + alt_seq[var.POS-start+1:]
            rsites = restriction_sites(ref_seq, alt_seq)
            if len(rsites) > 0:
                primer_set = primer3()
                for primer in primer_set.fetch_primers(ref_seq):
                    yield primer, restriction_sites, start, end
            

if __name__ == '__main__':
    args = docopt(__doc__, 
                  argv = debug,
                  options_first=False)

    v = vcf(args["<vcf>"])

    # Set VCF Options
    v.size = args["--size"]
    v.reference = args["--ref"]

    
    if args["--region"]:
        chrom, start, end = parse_region(args["--region"])
    else:
        chrom, start, end = [None]*3
    
    if args["--samples"]:
        if args["--samples"] == "ALL":
            samples = v.samples
        else:
            samples = args["--samples"].split(",")
            for sample in samples:
                if sample not in v.samples:
                    with indent(4):
                        puts_err(colored.red("\n" + sample + " Not found in VCF\n"))
                        exit()
    else:
        samples = None

    # Check for std. input
    if args["template"]:
        print('\t'.join(["CHROM_POS",
                         "REF",
                         "ALT",
                         "SAMPLE",
                         "EDIT_DISTANCE",
                         "SEQUENCE",
                         "0/0",
                         "0/1",
                         "1/1",
                         "./."]))
        for cvariant in v.fetch_variants_w_consensus(chrom, start, end, int(args["--size"]), args["--template"], samples):
            print(cvariant)


    elif args["indel"]:
        seq = v.fetch_variants_w_consensus(chrom, start, end)
        print(list(seq))
        exit()

    elif args["snipsnp"]:
        for primer, restriction_sites, start, end in v.extract_restriction():
            #print primer
            print restriction_sites
            print start, end

