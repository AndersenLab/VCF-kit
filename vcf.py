from subprocess import PIPE, Popen
from collections import OrderedDict
from tabulate import tabulate as tb
from ast import literal_eval
import re
import sys
import csv

_vcf_variable_types = {"Integer" : "Integer", "String": "Categorical", "Float" : "Float", "Flag": "Bool"}

standard_set_names = ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER"]
standard_set = [OrderedDict([("id","CHROM"), ("desc", "Chromosome"), ("type", "Categorical")]),
    OrderedDict([("id","POS"), ("desc", "Position"), ("type", "Integer")]),
    OrderedDict([("id","REF"), ("desc", "Reference Allele"), ("type", "Categorical")]),
    OrderedDict([("id","ALT"), ("desc", "Alternate Allele"), ("type", "Categorical")]),
    OrderedDict([("id","QUAL"), ("desc", "Variant Quality"), ("type", "Float")]),
    OrderedDict([("id","FILTER"), ("desc", "Filter"), ("type", "Categorical")])]


class bcolors:
    BOLD = "\033[1m"
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

def error(text):
    """ Reports an error to the user and exits """
    print("")
    print(bcolors.FAIL + "VCF-Toolbox Error " + bcolors.ENDC + text)
    print("")
    sys.exit(0)

def command(command):
  comm, err = Popen(command, stdout=PIPE, stderr=PIPE).communicate()
  if err != "":
    raise Exception(bcolors.WARNING + "BCFtools Error " + bcolors.ENDC + self.error)
  else:
    return comm.strip()    


def bc(text, color):
    return getattr(bcolors,color) + text + bcolors.ENDC




class vcf:
  def __init__(self, filename):
      # Start by storing basic information about the vcf/bcf and checking that an index exists.
      self.filename = filename
      self.header = command(["bcftools","view","-h",filename])
      
      # Samples
      self.samples = command(["bcftools","query","-l",filename]).split("\n")

      # Meta Data
      self.metadata = OrderedDict(re.compile(r'''^##(?P<key>[^<#]+?)=(?P<val>[^<#]+$)''', re.M).findall(self.header))

      # Contigs
      self.contigs = [x.split(",") for x in re.compile(r'''^##contig=<(?P<data>.*)>''', re.M).findall(self.header)]
      self.contigs = [{x.split("=")[0]:x.split("=")[1] for x in f} for f in self.contigs]
      
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
      self.format_set = {x["id"]:x for x in [m.groupdict() for m in r.finditer(self.header)]}
      self.info_format_variables = self.info_set.keys() + self.format_set.keys()
      self.complete_var_list = standard_set_names + self.info_format_variables

  def list_vars(self):
    # Resolve Variable Specification Issues

    print("")
    print(bc("STANDARD","BOLD"))
    print tb(standard_set, headers="keys", tablefmt="grid")
    print("")
    print(bc("INFO FIELDS","BOLD"))
    info = [OrderedDict([("id","INFO/" + x["id"]), ("desc", x["desc"]), ("type", x["type"])]) for x in self.info_set.values()]
    print tb(info, headers="keys", tablefmt="grid")
    print("")
    print(bc("FORMAT FIELDS","BOLD"))
    format = [OrderedDict([("id", "FORMAT/" + x["id"]), ("desc", x["desc"]), ("type", x["type"])]) for x in self.format_set.values()]
    print tb(format, headers="keys", tablefmt="grid")
    print("")
    print("""Variables have been prefixed with "INFO/" or "FORMAT/", \nhowever these prefixes are not necessary if the variable \nis unique within both groups.
        """)

  def resolve_variable(self, var):
    """ Finds variable in INFO or FORMAT Fields """
    if var.startswith("INFO/"):
        return "%" + var
    elif var in self.info_set.keys():
        return "%INFO/" + var
    else:
        return "[%" + var.replace("FORMAT/","") + "]"

  def format_var_name(self, var):
    """ Returns proper specification of variable """
    var = var.replace("%", "")
    search_var = var.replace("FORMAT/", "").replace("INFO/","")
    if search_var not in self.complete_var_list:
        error("%s not found; Use listvars to see available variables" % var)
    elif self.info_format_variables.count(var) > 1:
        error("%s is defined in FORMAT and INFO sets; prefix with INFO/ or FORMAT/" % var)
    elif var.startswith("INFO/") and var.replace("INFO/", "") not in self.info_set.keys():
        error("%s not found in INFO variables" % var)
    elif var in standard_set_names:
        return "%" + var
    else:
        return self.resolve_variable(var)
  
  def get_variable_type(self, var):
    pass

  def query(self, x, y=None, region=None, include=None):
    x = self.format_var_name(x)
    q = ["bcftools","query","-f",'\"%s\\n\"' % x, self.filename]
    print q
    print " ".join(q)
    

  def _parse_stats(lines):
    stats = {}
    values_lines = []
    lines = lines.split("\n")
    for l in lines:
      if l.startswith("#"):
        if len(values_lines) > 0:
          stats[stat_section] = values_lines
          values_lines = []
        parse_header = re.sub("\[[0-9]+\]","",l.replace("# ","")).split("\t")
        stat_section = parse_header[0]
        columns = parse_header[1:]
      else:
        formatted_values = coerce_types(l.split("\t"))[1:]
        values_lines.append(zip(columns,formatted_values))
    return stats
