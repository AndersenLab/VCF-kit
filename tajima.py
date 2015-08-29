#! /usr/bin/env python
"""
usage:
  tb.py tajima 

Example

command:
  tajima        Calculate Tajima's D

options:
  -h --help                   Show this screen.
  --version                   Show version.



"""
from docopt import docopt
from subprocess import call
from utils.vcf import *
import sys



debug = None
if len(sys.argv) == 1:
    debug = ['tajima']


class tajima(vcf):
  """
    Subclass of the vcf object
    used for calculating tajima's D
  """
  def __init__(self, filename):
    vcf.__init__(self,filename)

    # Tajima D Constants
    self.a1 = sum([1.0/d for d in xrange(1,self.n-1)])
    self.a2 = sum([1.0/d**2 for d in xrange(1,self.n-1)])
    self.b1 = (self.n + 1.0) / (3.0*(self.n - 1))
    self.b2 = ( 2 * (self.n**2 + self.n + 3.0) / 
              (9.0 * self.n * (self.n -1.0)))
    self.c1 = self.b1 - (1.0 / self.a1)
    self.c2 = self.b2 - (self.n + 2 / (self.a1*self.n)) + (self.a2 / (self.a1**2) )
    self.e1 = self.c1 / self.a1
    self.e2 = self.c2 / (self.a1**2 + self.a2)

  def calc_tajima(self):
      for variant_interval in self.window(window_size=1000000, shift_method="POS-Interval"):
          self.gt = np.vstack([x.gt_types for x in variant_interval])
          self.segregating_sites = sum([np.any(p) for p in np.equal(3, self.gt)])
          print(self.segregating_sites)




class variant_ops:
    """
        Variant operator - takes a variant interval as input
        and enables certain operations to be performed
    """
    def __init__(self, variant_interval, vcf):
        # self.shift_method = shift_method
        self.interval = variant_interval.interval()
        self.gt = np.vstack([x.gt_types for x in variant_interval])
        self.n = vcf.n 
        self.segregating_sites = sum([np.any(p) for p in np.equal(3, self.gt)])
        #print self.gt.T
        #print np.vstack([x.gt_bases for x in variant_interval]).T
        print self.segregating_sites
        self.a1 = H(self.n) # Number of DNA sequences



if __name__ == '__main__':
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug,
                  options_first=True)
    print(args)
    print "TAJIMA"


v = tajima("test.vcf.gz")
print dir(v)
v.calc_tajima()


