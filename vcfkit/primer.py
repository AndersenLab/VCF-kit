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
  --size=<int>                Amplicon size/2 (Upstream and downstream length) [default: 50]
  --box-variants              Add second column for the sequence with the variant boxed.


"""
from docopt import docopt
from utils import message
from utils.primer_vcf import *
from utils.reference import *
from utils.fasta import *
import sys
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


if __name__ == '__main__':
    args = docopt(__doc__,
                  argv=debug,
                  options_first=False)

    # Ensure user has specified a reference.
    if args["--ref"] is None:
        exit(message("Must specify a reference with --ref", color="red"))

    v = primer_vcf(args["<vcf>"], reference=args["--ref"])

    # Region
    if args["--region"]:
        v.region = args["--region"]
    else:
        v.region = None

    # Samples
    if args["--samples"]:
        if args["--samples"] == "ALL":
            v.output_samples = v.samples
        else:
            v.output_samples = args["--samples"].split(",")
            for sample in v.output_samples:
                if sample not in v.samples + ["REF", "ALT", ]:
                    exit(message(sample + " not found in VCF", "red"))

    v.box_variants = args["--box-variants"]
    v.size = int(args["--size"])
    v.template = args["--template"]

    # Check for std. input
    if args["template"]:
        v.mode = "template"
        v.print_primer_header()
        for cvariant in v.fetch_template():
            print(cvariant)

    elif args["indel"]:
        v.mode = "indel"
        for cvariant in v.fetch_indel_primers():
            print(list(cvariant))


    elif args["snpsnp"]:
        for primer, restriction_sites, start, end in v.extract_restriction():
            # print primer
            print restriction_sites
            print start, end
