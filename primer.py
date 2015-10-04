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


    def extract_restriction(self, window = 20000):
        for varset in self.window(window_size=window, step_size = 1, shift_method="POS-Sliding"):
            print len(varset), varset
            start = varset.lower_bound
            end = varset.upper_bound
            print start, end
            print varset[0].ALT, varset[0].POS
            CHROM = varset[0].CHROM
            ref = self.reference[CHROM][start:end].seq
            alt = ref
            for var in varset:
                print var.CHROM, var.POS
                # Spike alt sequence with all variants.
                alt = alt[:var.POS-start] + var.ALT[0] + alt[var.POS-start+1:]
            ref = Seq(ref, DNA_SET)
            alt = Seq(alt, DNA_SET)
            ref = dict(AllEnzymes.search(ref).items())
            alt = dict(AllEnzymes.search(alt).items())
            ref = {k:v for k,v in ref.items() if len(v) > 0 and len(v) <= 3 and
                                                             ref[k] != alt[k] and
                                                             abs(len(ref[k]) - len(alt[k])) == 1}
            alt = {k:v for k,v in alt.items() if k in ref.keys()}
            yield varset, ref, alt, start, end
            


if __name__ == '__main__':
    #print debug
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug,
                  options_first=False)
    print args
    # Locate Reference
    v = seq_vcf(args["<vcf>"], args["--ref"])
    print dir(seq_vcf)
    if args["snpsnp"] == True:
        for varset, ref, alt, start, end in v.extract_restriction():
            print varset, ref, alt, start, end



    x = seq_vcf(args["<vcf>"])
