#! /usr/bin/env python
"""
usage:
  tb.py primer [--ref=<reference>] <vcf>

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


    def extract(self):
        for i in self:
            window = 1000
            start = i.POS - window + 1
            end = i.POS + window
            yield self.reference[i.CHROM][start:end]
            


if __name__ == '__main__':
    #print debug
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug,
                  options_first=False)
    
    # Locate Reference
    v = seq_vcf(args["<vcf>"], args["--ref"])
    for i in v.extract():
        print i
        exit()

    search_results2 = AllEnzymes.search( Seq("tattgaaaaaaac", DNA_SET ))
    search_results = AllEnzymes.search(  Seq("tattgaagtaaac", DNA_SET ))
    sr2 = search_results2.items()
    for i in search_results2:
      print i.size, i.elucidate()
    #print search_results2.items()
    print dir(AllEnzymes)
    print [(k,v) for k,v in search_results.iteritems() if len(v) > 0 and (k,v) not in sr2]
    

    x = seq_vcf(args["<vcf>"])
