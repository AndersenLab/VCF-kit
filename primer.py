#! /usr/bin/env python
"""
usage:
  tb.py primer [--ref=<reference>] <vcf>
  tb.py primer snpsnp [--ref=<reference>] <vcf>

Example

options:
  -h --help                   Show this screen.
  --version                   Show version.



"""
from docopt import docopt
from clint.textui import colored, puts, indent
from subprocess import call
from utils.vcf import *
from utils.fasta import *
import sys
import gc
import os
from glob import glob
from Bio.Seq import Seq 
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as DNA_SET
from Bio.Restriction import AllEnzymes
import Bio

debug = None
if len(sys.argv) == 1:
    debug = ['primer', "--ref=WBcel235", "test.vcf.gz"]


class seq_vcf(vcf):
    """
        Subclass of vcf that extracts genomic context of variants from reference.
    """
    def __init__(self, filename, reference):
        vcf.__init__(self,filename)
        if not os.path.isfile(reference):
            module_path = os.path.split(os.path.realpath(__file__))[0]
            self.reference = module_path + "/reference/" + args["--ref"] + ".fa.gz"
            if not os.path.isfile(self.reference):
                with indent(4):
                    puts(colored.red("\nGenome does not exist!\n"))
                    puts(colored.blue("Available Genomes:\n"))
                with indent(8):
                    for i in glob(module_path + "/reference/*.fa.gz"):
                        genome = os.path.split(i)[1].replace(".fa.gz", "")
                        puts(genome)
                with indent(4):
                    puts(colored.blue("\nDownload more with tb.py genome\n"))
                    exit()
            self.reference = Fasta(self.reference)
            if len([x for x in self.contigs.keys() if x not in self.reference.keys()]) > 0:
                with indent(4):
                    puts(colored.yellow("""\nWarning: VCF / Reference contigs don't match\nUsing contig names within VCF.\n"""))
                    self.reference.alt_contig_names = dict(zip(self.contigs.keys(), self.reference.keys(),))


    def extract_ref_alt(self):
        for i in self:
            if i.is_snp:
                window = 1000
                start = i.POS - window + 1
                end = i.POS + window
                ref = self.reference[i.CHROM][start:i.POS-1].seq + i.REF + self.reference[i.CHROM][i.POS+1:end].seq
                alt = self.reference[i.CHROM][start:i.POS-1].seq + i.ALT[0] + self.reference[i.CHROM][i.POS+1:end].seq
                yield Seq(ref, DNA_SET), Seq(alt, DNA_SET)
            


if __name__ == '__main__':
    #print debug
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug,
                  options_first=False)
    print args
    # Locate Reference
    v = seq_vcf(args["<vcf>"], args["--ref"])

    if args["snpsnp"] == True:
        for ref, alt in v.extract_ref_alt():
            restriction_ref = dict(AllEnzymes.search(ref).items())
            restriction_alt = dict(AllEnzymes.search(alt).items())
            ref_sites = {k:v for k,v in restriction_ref.items() if len(v) > 0 and len(v) <= 3 and
                                                                 restriction_ref[k] != restriction_alt[k] and
                                                                 abs(len(restriction_ref[k]) - len(restriction_alt[k])) == 1}
            alt_sites = {k:v for k,v in restriction_alt.items() if k in ref_sites.keys()}
            print ref_sites, alt_sites



    x = seq_vcf(args["<vcf>"])
