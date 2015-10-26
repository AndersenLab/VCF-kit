from subprocess import Popen, PIPE
from utils import *
from collections import defaultdict
from pprint import pprint as pp
import csv


field = ["query_id",
         "CHROM",
         "perc_identity",
         "alignment length",
         "mismatches",
         "gap opens",
         "q_start",
         "q_end",
         "s_start",
         "s_end",
         "evalue",
         "bit_score"]



class blastn:

    def __init__(self, db):
        self.db = db

    def blast(self, q, chrom = None, pos = None):
        blastn_query = "echo {q} | blastn -query - -db={self.db} -outfmt=6 -evalue 1 -word_size=7".format(**locals())
        self.query = q
        self.query_length = len(q)
        self.chrom = chrom
        self.pos = pos
        resp, err = Popen(blastn_query,
                          stdout=PIPE,
                          stderr=PIPE,
                          shell=True).communicate()

        if err != "":
            raise Exception(err)


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


#blaster = blastn("reference/WS245/WS245.fa")
#print pp(blaster.blast("ttaggcttaggcttaggc", chrom = "X", pos = 22))
#print blaster.query_stat()
