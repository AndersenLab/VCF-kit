from subprocess import Popen, PIPE
from vcfkit.utils import *
from collections import defaultdict, OrderedDict
from pprint import pprint as pp
import csv
from clint.textui import puts_err
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from copy import copy
from cStringIO import StringIO
from vcfkit.utils import *

def clamp(n, minn, maxn):
    return max(min(maxn, n), minn)

def fastq_mean(x):
    if len(x) > 0:
        return round(sum(x)*1.0/len(x), 2)
    else:
        return 0


class blast_variant:
    """
        Class representing a blast variant
    """
    output_order = ["CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "seq_gt",
                    "vcf_gt",
                    "sample",
                    "variant_type",
                    "classification",
                    "index",
                    "alignment_start",
                    "alignment_end",
                    "strand",
                    "context",
                    "gaps",
                    "mismatch",
                    "evalue",
                    "bitscore",
                    "phred_quality",
                    "phred_quality_window",
                    "description"]

    def __init__(self,
                 blast_result,
                 POS,
                 ref_out,
                 alt_out,
                 index,
                 context,
                 gaps,
                 mismatch,
                 phred_quality = "",
                 phred_quality_window = ""):

        self.CHROM = blast_result["sacc"]
        self.POS = POS
        posp1 = self.POS + 1
        self.CHROM_POS = "{self.CHROM}:{self.POS}-{posp1}".format(**locals())
        self.classification = None
        self.REF = ref_out.strip("-")
        self.ALT = None
        self.seq_gt = alt_out.strip("-")
        self.is_variant = (self.REF != self.seq_gt)
        self.vcf_gt = ""
        self.alignment_start = blast_result["sstart"]
        self.alignment_end = blast_result["send"]
        self.strand = blast_result["strand"]
        self.context = context
        self.gaps = gaps
        self.mismatch = mismatch
        self.index = index
        self.evalue = blast_result["evalue"]
        self.bitscore = blast_result["bitscore"]
        if self.REF != self.seq_gt:
            if len(self.REF) == len(self.seq_gt):
                self.variant_type = "snp"
            elif len(self.REF) > len(self.seq_gt):
                self.variant_type = "deletion"
            elif len(self.REF) < len(self.seq_gt):
                self.variant_type = "insertion"
        else:
            self.variant_type = ""
        self.phred_quality = phred_quality
        self.phred_quality_window = phred_quality_window

    def classify_variant(self):
        # Determines whether variant is TP, TN, FP, FN
        pass

    def chrom_pos_allele(self):
        """
            Simplified representation of variant
        """
        return [self.CHROM, self.POS, self.seq_gt]

    def region(self):
        return "{self.CHROM}:{self.alignment_start}-{self.alignment_end}".format(**locals())

    def __str__(self):
        return '\t'.join([str(getattr(self, x)) for x in blast_variant.output_order])


class blast:
    def __init__(self, db, num_alignments=1, word_size = 28):
        self.db = db
        self.num_alignments = num_alignments
        self.word_size = word_size
        self.output_format = ["sacc",
                              "pident",
                              "gaps",
                              "mismatch",
                              "length",
                              "qstart",
                              "qend",
                              "sstart",
                              "send",
                              "evalue",
                              "bitscore",
                              "qseq",
                              "sseq"]
        self.output_string = ' '.join(self.output_format)
        self.blastn_query_str = ' '.join(["echo {self.query} | blastn -query - -db={self.db} ",
                             "-outfmt '6 {self.output_string}'",
                             "-num_alignments {self.num_alignments}",
                             "-word_size {self.word_size}"])

    def print_alignment(self, q, start=0, end=None):
        """ 
            Print the alignment with matches

            * [ ] Add qstart/qend; sstart/ssend
        """
        blast_result = self.blast_search(q)
        if end is None:
            end = len(blast_result["sseq"])
        sseq = blast_result["sseq"]
        qseq = blast_result["qseq"]
        spacer = []
        for i in xrange(len(sseq)):
            if qseq[i] == sseq[i]:
                spacer.append("|")
            else:
                spacer.append(" ")
        spacer = ''.join(spacer)
        print(blast_result["sseq"][start:end])
        print(spacer[start:end])
        print(blast_result["qseq"][start:end])

    def blast_search(self, q):
        """
            Perform a blast search and store the result
        """
        if type(q) == SeqRecord:
            self.query = q.seq
            if hasattr(q, "letter_annotations"):
                if 'phred_quality' in q.letter_annotations.keys():
                    self.query_qual = q.letter_annotations["phred_quality"]
                else:
                    self.query_qual = None
        else:
            self.query = q

        self.query_length = len(q)
        blastn_query = self.blastn_query_str.format(**locals())
        resp, err = Popen(blastn_query,
                          stdout=PIPE,
                          stderr=PIPE,
                          shell=True).communicate()
        if not resp:
           return None
        if err:
           raise Exception(err)

        # Format variables
        resp = [OrderedDict(zip(self.output_format,
                                map(autoconvert, x.split("\t")))) 
                for x in resp.splitlines()]
        return resp

    def check_primer(self, q):
        blast_results = self.blast_search(q)
        return len([i for i in blast_results if i['length'] > 14])


    def blast_call(self, q):
        """
            Call variants from blast alignments
        """
        blast_result = self.blast_search(q)
        if blast_result is None:
            yield None
        else:
            CHROM = blast_result["sacc"]
            start = blast_result["sstart"]
            end = blast_result["send"]
            ref = blast_result["sseq"]
            alt = blast_result["qseq"]
            gaps = blast_result["gaps"]
            mismatch = blast_result["mismatch"]

            # Compare fasta sites
            ref_out = ""
            alt_out = ""
            len_insertions = 0
            len_deletions = 0
            phred_quality = ""
            phred_quality_window = ""
            i = 0
            while i < len(ref)-1:
                if alt[i+1] == "-":
                    ref_out = ref[i]
                    alt_out = alt[i]
                    while alt[i+1] == "-":
                        i += 1
                        ref_out += ref[i]
                        alt_out += alt[i]
                    len_deletions += len(alt_out) - 1
                elif ref[i+1] == "-":
                    ref_out = ref[i]
                    alt_out = alt[i]
                    while ref[i+1] == "-":
                        i += 1
                        ref_out += ref[i]
                        alt_out += alt[i]
                    len_insertions += len(ref_out) - 1
                else:
                    ref_out, alt_out = ref[i], alt[i]
                

                POS = i - len_insertions
                context_start = clamp(i-10,0,len(ref))
                context_end = clamp(i+10,0,len(ref))
                context = alt[context_start:i] + "[" + alt_out + "]" + alt[i+1:context_end]
                context = context.replace("-","")
                index = i + blast_result["qstart"] - len_deletions - len(alt_out)

                if 'qqual' in blast_result:
                    window_start = clamp(i-6, 0, len(blast_result["qqual"]))
                    window_end = clamp(i+7, 0, len(blast_result["qqual"]))
                    phred_quality_window = fastq_mean([x for x in blast_result["qqual"][window_start:window_end] if x != "NA"])
                    phred_quality = blast_result["qqual"][i]

                variant_out = blast_variant(blast_result=blast_result,
                                            POS=POS + start,
                                            ref_out=ref_out,
                                            alt_out=alt_out,
                                            index=index,
                                            context=context,
                                            gaps=gaps,
                                            mismatch=mismatch,
                                            phred_quality=phred_quality,
                                            phred_quality_window=phred_quality_window)
                yield variant_out
                i += 1

