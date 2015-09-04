from subprocess import check_output
from collections import OrderedDict

class Fasta:
    """
        Extract regions from fasta files indexed with samtools
    """
    def __init__(self, reference):
        self.reference = reference

    def __getitem__(self, chrom_pos):
        print chrom_pos, "CHROM"
        if isinstance(chrom_pos, slice):
            start = chrom_pos.start
            if start < 1:
                start = 1
            end = chrom_pos.stop
            chrom_pos = "{self.chrom}:{start}-{end}".format(**locals())
            query_string = ["samtools", "faidx", self.reference, chrom_pos]
            return check_output(query_string)
        self.chrom = chrom_pos
        return self

    def keys(self, weight = False):
        with open(self.reference + ".fai") as f:
            chrom_length = [x.strip().split("\t")[0:2] for x in f.readlines()]
            chrom_length = [[x[0],int(x[1])] for x in chrom_length]
            if weight == True:
                genome_length = sum([x[1] for x in chrom_length])
                return OrderedDict([x[0],1.0*x[1]/genome_length] for x in chrom_length)
            else:
                return OrderedDict(chrom_length)
