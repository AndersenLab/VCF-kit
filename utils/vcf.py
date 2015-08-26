from subprocess import PIPE, Popen
from collections import OrderedDict, deque
import re



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


class vcf:
    def __init__(self, filename):
        self.filename = filename
        self.__reader = Popen(["bcftools", "view", "-r", "I:1-10000", self.filename],
                                         stdout=PIPE,
                                         stderr=PIPE).stdout
        self.header = self.__read_header()
        self.samples = self.header.strip().split("\n")[-1].split("\t")[9:]
        # Meta Data
        self.metadata = OrderedDict(re.compile(r'''^##(?P<key>[^<#]+?)=(?P<val>[^<#]+$)''', re.M).findall(self.header))

        # Contigs
        self.contigs = [x.split(",") for x in re.compile(r'''^##contig=<(?P<data>.*)>''', re.M).findall(self.header)]
        self.contigs = [dict([(x.split("=")[0],autoconvert(x.split("=")[1])) for x in f]) for f in self.contigs]
        self.contigs = dict([(x["ID"],x) for x in self.contigs])
        
        # Info
        r = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.info_set = {x["id"]:x for x in [m.groupdict() for m in r.finditer(self.header)]}
        
        # Filter
        r = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.filter_set = {x["id"]:x for x in [m.groupdict() for m in r.finditer(self.header)]}
        
        # Format
        r = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)

    def __read_header(self):
        header_lines = []
        while True:
            line = self.__reader.next().strip()
            if line.startswith("#"):
                header_lines.append(line)
            elif line.startswith("##"):
                header_lines.append(line)
            else:
                self.__first_line = line
                return "\n".join(header_lines)

    def next(self):
        if self.__first_line is not None:
            first_line = self.__first_line
            self.__first_line = None
            return first_line
        else:
            return self.__reader.next().strip()

    def window(self, windowsize, shift_type, shift_kind):
        # Add start and end
        line = "n"
        result_list = deque(maxlen = windowsize)
        try:
            while True:
                if shift_kind == "SNP":
                    if shift_type == "sliding":
                        result_list.append(self.next())
                    elif shift_type == "interval":
                        result_list = []
                        for snp in xrange(0,windowsize):
                            result_list.append(self.next())
                elif shift_kind == "position":
                    pass
                result = [int(x.split("\t")[1]) for x in result_list]
                if len(result) == windowsize:
                    yield result 
        except StopIteration:
            pass
                






