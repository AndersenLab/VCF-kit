#! /usr/bin/env python
"""
usage:
  tb.py primer --option <vcf>

Example

options:
  -h --help                   Show this screen.
  --version                   Show version.



"""
from docopt import docopt
from subprocess import call
from utils.vcf import *
import sys
import gc

from Bio.Seq import Seq 
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as DNA_SET
from Bio.Restriction import AllEnzymes
import Bio

debug = None
if len(sys.argv) == 1:
    debug = ['primer', "--option", "test.vcf.gz"]



if __name__ == '__main__':
    print debug
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug,
                  options_first=False)
    v = vcf(args["<vcf>"])

    search_results2 = AllEnzymes.search( Seq("tattgaaaaaaac", DNA_SET ))
    search_results = AllEnzymes.search(  Seq("tattgaagtaaac", DNA_SET ))
    sr2 = search_results2.items()
    for i in search_results2:
      print i.size, i.elucidate()
    #print search_results2.items()
    print dir(AllEnzymes)
    print [(k,v) for k,v in search_results.iteritems() if len(v) > 0 and (k,v) not in sr2]
    

    x = vcf(args["<vcf>"])
