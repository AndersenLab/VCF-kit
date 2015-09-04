from subprocess import check_output

class Fasta:
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

    def keys(self):
        with open(self.reference + ".fai") as f:
            return [x.split("\t")[0] for x in f.read().strip().split("\n")]
