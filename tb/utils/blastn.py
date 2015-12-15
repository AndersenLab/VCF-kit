from subprocess import Popen, PIPE
from utils import *
from collections import defaultdict, OrderedDict
from pprint import pprint as pp
import csv
from clint.textui import puts_err
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from cStringIO import StringIO


def boolify(s):
    if s == 'True':
        return True
    if s == 'False':
        return False
    raise ValueError("huh?")

def autoconvert(s):
    for fn in (boolify, int, float):
        try:
            return fn(s)
        except ValueError:
            pass
    return s

def clamp(n, minn, maxn):
    return max(min(maxn, n), minn)

def fastq_mean(x):
    if len(x) > 0:
        return round(sum(x)*1.0/len(x), 2)
    else:
        return 0

class blast:
    def __init__(self, db, num_alignments=1):
        self.db = db
        self.num_alignments = num_alignments
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
                             "-num_alignments {self.num_alignments}"])

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
        self.blastn_query = self.blastn_query_str.format(**locals())
        resp, err = Popen(self.blastn_query,
                          stdout=PIPE,
                          stderr=PIPE,
                          shell=True).communicate()

        if not resp:
           return None
        if err:
           raise Exception(err)

        # Format variables
        resp = map(autoconvert, [x.strip() for x in resp.split("\n")[0].split("\t")])
        # Convert to dictionary
        resp = OrderedDict(zip(self.output_format, resp))
        # Append quality values if available, and insert gaps
        if self.query_qual:
            qseq = resp["qseq"]
            query_qualities = self.query_qual[resp["qstart"] - 1:resp["qend"]]
            # Add gaps
            inserted_gaps = 0
            for i in xrange(len(qseq)):
                if qseq[i] == "-":
                    query_qualities.insert(i - inserted_gaps , "NA")
                    inserted_gaps += 1
            resp["qqual"] = query_qualities
        return resp


    def call_header(self, qual = False, sample = False):
        """
            Print output header for blast_call
        """
        self.header = ["CHROM",
                       "POS",
                       "ref",
                       "gt",
                       "variant_type",
                       "alignment_start",
                       "alignment_end",
                       "gaps",
                       "mismatch",
                       "index",
                       "context"]
        if sample:
            self.header += ["sample", "description"]
        if qual:
            self.header.insert(7, "quality_window")
            self.header.insert(7, "quality")
        print '\t'.join(self.header)


    def blast_call(self, q, all_sites = False):
        """
            Call variants from blast alignments
        """
        r = self.blast_search(q)
        if r is None:
            yield None
        else:
            chrom = r["sacc"]
            start = r["sstart"]
            end = r["send"]
            ref = r["sseq"]
            alt = r["qseq"]
            gaps = r["gaps"]
            mismatch = r["mismatch"]

            # Compare fasta sites
            ref_out = ""
            alt_out = ""
            len_insertions = 0
            len_deletions = 0
            i = 0
            while i < len(ref)-1:
                if alt[i+1] == "-":
                    variant_type = "deletion"
                    ref_out = ref[i]
                    alt_out = alt[i]
                    while alt[i+1] == "-":
                        i += 1
                        ref_out += ref[i]
                        alt_out += alt[i]
                    len_deletions += len(alt_out) - 1
                elif ref[i+1] == "-":
                    variant_type = "insertion"
                    ref_out = ref[i]
                    alt_out = alt[i]
                    while ref[i+1] == "-":
                        i += 1
                        ref_out += ref[i]
                        alt_out += alt[i]
                    len_insertions += len(ref_out) - 1
                else:
                    ref_out, alt_out = ref[i], alt[i]
                    if ref_out != alt_out:
                        variant_type = "snp"
                    else:
                        variant_type = ""
                if ref_out != alt_out or (ref_out == alt_out and all_sites):
                    POS = i - len_insertions
                    context_start = clamp(i-10,0,len(ref))
                    context_end = clamp(i+10,0,len(ref))
                    context = alt[context_start:i] + "[" + alt[i] + "]" + alt[i+1:context_end]
                    context = context.replace("-","")
                    out_result = [chrom,
                                  POS + start,
                                  ref_out.strip("-"),
                                  alt_out.strip("-"),
                                  variant_type,
                                  start,
                                  end,
                                  gaps,
                                  mismatch,
                                  i + r["qstart"] - len_deletions - len(alt_out),
                                  context]
                    if 'qqual' in r:
                        window_start = clamp(i-6, 0, len(r["qqual"]))
                        window_end = clamp(i+7, 0, len(r["qqual"]))
                        window = [x for x in r["qqual"][window_start:window_end] if x != "NA"]
                        out_result.insert(7, str(fastq_mean(window)))
                        out_result.insert(7, str(r["qqual"][i]))
                    yield map(str,out_result)
                i += 1




"""
OLD - remove/rewrite for snip-snp

        resp = [dict(zip(field,map(autoconvert,x.split("\t")))) for x in resp.strip().split("\n")]
        exact_matches = 0
        for i in resp:
            i["query_length"] = abs(i["q_end"] - i["q_start"]) + 1
            i["subject_length"] = abs(i["s_end"] - i["s_start"]) + 1
            i["query_string"] = q[i["q_start"]-1:i["q_end"]]

            # Determine if its a
            i["perfect_match"] = (i["query_string"] == self.query)
            i["bp_mismatch"] = (i["subject_length"] == self.query_length
                                and i["mismatches"] == 1)
            i["end_mismatch"] = (self.query.find(i["query_string"]) >= 0
                                 and i["perfect_match"] == False)
            if chrom is not None:
                i["perfect_match_same_chrom"] = (i["perfect_match"] is True
                                                 and i["CHROM"] == self.chrom)
            else:
                i["perfect_match_same_chrom"] = False
            if chrom is not None and pos is not None:
                i["nearby_matches"] = ((i["perfect_match"] is True or
                                        i["bp_mismatch"] is True or
                                        i["end_mismatch"] is True) and
                                       i["CHROM"] == self.chrom and
                                       abs(i["s_start"] - pos) < 10000)
            else:
                i["nearby_matches"] = False

        self.resp = resp

        return self

    def query_stat(self):
        self.stats = {}
        for i in ["perfect_match", "bp_mismatch", "end_mismatch", "perfect_match_same_chrom", "nearby_matches"]:
            self.stats[i] = len([x for x in self.resp if x[i] is True])
        return self.stats
"""
