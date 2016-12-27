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
import sys
import signal
signal.signal(signal.SIGINT, lambda x,y: sys.exit(0))



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
                    "reference",
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
                 reference_seq,
                 query_seq,
                 seq_out,
                 index,
                 reference_index,
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
        self.reference = reference_seq[reference_index]
        self.seq_gt = seq_out.strip("-")
        self.vcf_gt = ""
        self.alignment_start = blast_result["sstart"]
        self.alignment_end = blast_result["send"]
        self.strand = {"plus": "+", "minus": "-"}[blast_result["sstrand"]]
        self.context = context
        self.gaps = gaps
        self.mismatch = mismatch
        self.index = index
        self.evalue = blast_result["evalue"]
        self.bitscore = blast_result["bitscore"]
        self.phred_quality = phred_quality
        self.phred_quality_window = phred_quality_window

    def fetch_variant_type(self):
        if self.REF != self.ALT:
            if len(self.REF) == len(self.seq_gt):
                variant_type = "snp"
            elif len(self.REF) > len(self.seq_gt):
                variant_type = "deletion"
            elif len(self.REF) < len(self.seq_gt):
                variant_type = "insertion"
        else:
            variant_type = ""
        self.variant_type = variant_type

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
                              "sseq",
                              "sstrand"]
        self.output_string = ' '.join(self.output_format)
        self.blastn_query_str = ' '.join(["echo {self.query} | blastn -query - -db={self.db} ",
                             "-outfmt '6 {self.output_string}'",
                             "-num_alignments {self.num_alignments}",
                             "-word_size {self.word_size}"])

    def fetch_reference_seq(self, chrom, start, end):
        end += 20
        fetch_ref_cmd = ["samtools",
         "faidx",
         self.db,
         "{chrom}:{start}-{end}".format(**locals())]
        resp, err = Popen(fetch_ref_cmd,
               stdout=PIPE,
               stderr=PIPE).communicate()
        return ''.join(resp.splitlines()[1:]).upper()


    def print_alignment(self, q, start=0, end=None):
        """ 
            Print the alignment with matches
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
        if blast_results:
            return len([i for i in blast_results if i['length'] > 14])
        else:
            return 0


    def blast_call(self, q):
        """
            Call variants from blast alignments
        """
        blast_result = self.blast_search(q)
        if blast_result is None:
            yield None
        else:
            for bresult in blast_result:
                CHROM = bresult["sacc"]
                start = bresult["sstart"]
                end = bresult["send"]
                ref = bresult["sseq"]
                alt = bresult["qseq"]
                gaps = bresult["gaps"]
                mismatch = bresult["mismatch"]

                # Retrieve reference sequence
                reference_seq = self.fetch_reference_seq(CHROM, start, end)

                # Compare fasta sites
                ref_out = ""
                seq_out = ""
                len_insertions = 0
                len_deletions = 0
                phred_quality = ""
                phred_quality_window = ""
                i = 0
                while i < len(ref)-1:
                    if alt[i+1] == "-":
                        ref_out = ref[i]
                        seq_out = alt[i]
                        while alt[i+1] == "-":
                            i += 1
                            ref_out += ref[i]
                            seq_out += alt[i]
                        len_deletions += len(seq_out) - 1
                    elif ref[i+1] == "-":
                        ref_out = ref[i]
                        seq_out = alt[i]
                        while ref[i+1] == "-":
                            i += 1
                            ref_out += ref[i]
                            seq_out += alt[i]
                        len_insertions += len(ref_out) - 1
                    else:
                        ref_out, seq_out = ref[i], alt[i]
                    

                    POS = i - len_insertions
                    context_start = clamp(i-10,0,len(ref))
                    context_end = clamp(i+10,0,len(ref))
                    context = alt[context_start:i] + "[" + seq_out + "]" + alt[i+1:context_end]
                    context = context.replace("-","")
                    index = i + bresult["qstart"] - len_deletions - len(seq_out)

                    if 'qqual' in bresult:
                        window_start = clamp(i-6, 0, len(bresult["qqual"]))
                        window_end = clamp(i+7, 0, len(bresult["qqual"]))
                        phred_quality_window = fastq_mean([x for x in bresult["qqual"][window_start:window_end] if x != "NA"])
                        phred_quality = bresult["qqual"][i]

                    variant_out = blast_variant(blast_result=bresult,
                                                POS=POS + start,
                                                reference_index=POS,
                                                index=index,
                                                seq_out=seq_out,
                                                query_seq=q,
                                                reference_seq=reference_seq,
                                                context=context,
                                                gaps=gaps,
                                                mismatch=mismatch,
                                                phred_quality=phred_quality,
                                                phred_quality_window=phred_quality_window)
                    yield variant_out
                    i += 1

