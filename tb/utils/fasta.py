from subprocess import check_output
from collections import OrderedDict


class sequence:
    """
        Sequence object for keeping track
        of chromosomal position.
    """
    def __init__(self, chrom, start, end, seq):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.seq = seq

    def __repr__(self):
        seq_repr = ">{self.chrom}:{self.start}-{self.end}\n{self.seq}"
        return seq_repr.format(**locals())


class Fasta:
    """
        Extract regions from fasta files indexed with samtools
    """
    def __init__(self, reference):
        self.reference = reference
        self.alt_contig_names = {}

    def __getitem__(self, chrom_pos):
        if isinstance(chrom_pos, slice):
            start = chrom_pos.start
            if start < 1:
                start = 1
            end = chrom_pos.stop
            chrom_pos = "{self.chrom}:{start}-{end}".format(**locals())
            query_string = ["samtools", "faidx", self.reference, chrom_pos]
            seq = ''.join(check_output(query_string).strip().split("\n")[1:])
            seq = sequence(self.chrom_name, start, end, seq)
            return seq
        self.chrom = chrom_pos
        self.chrom_name = chrom_pos
        if chrom_pos in self.alt_contig_names.keys():
            self.chrom = self.alt_contig_names[chrom_pos]
        return self

    def keys(self, weight=False):
        with open(self.reference + ".fai") as f:
            chrom_length = [x.strip().split("\t")[0:2] for x in f.readlines()]
            chrom_length = [[x[0], int(x[1])] for x in chrom_length]
            if weight is True:
                genome_length = sum([x[1] for x in chrom_length])
                return OrderedDict([x[0], 1.0*x[1]/genome_length]
                                   for x in chrom_length)
            else:
                return OrderedDict(chrom_length)
