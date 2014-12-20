from subprocess import PIPE, Popen, call
from collections import OrderedDict
from tabulate import tabulate as tb
from ast import literal_eval
from datetime import datetime
from utils import *
import shlex
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
standard_set_types = {"CHROM":"Categorical",
                      "POS":"Integer",
                      "REF":"Categorical",
                      "ALT":"Categorical",
                      "QUAL":"Float",
                      "FILTER":"Categorical"}



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


  def format_var_for_query(self, var):
    """ Returns proper specification of variable,
    variable name for plotting, and split variables
    as needed.
    """
    if var == None:
        return {"query": None, "df": None , "type": None}
    var = var.replace("%", "")
    search_var = var.replace("FORMAT/", "").replace("INFO/","")
    if search_var not in self.complete_var_list:
        error("%s not found; Use listvars to see available variables" % var)
    elif self.info_format_variables.count(var) > 1:
        error("%s is defined in FORMAT and INFO sets; prefix with INFO/ or FORMAT/" % var)
    elif var.startswith("INFO/") and var.replace("INFO/", "") not in self.info_set.keys():
        error("%s not found in INFO variables" % var)
    elif var in standard_set_names:
        return {"query":"%" + var, "df": var, "type": standard_set_types[var], "number": 1}
    else:
        if var.startswith("INFO/") or var in self.info_set.keys():
            var_info = self.info_set[var.replace("INFO/","")]
            var_number = int(var_info["number"])
            if var_number > 1:
                var_data_rep =  ','.join([var + "." + str(v) for v in range(0,int(var_number))])
            else:
                var_data_rep = var
            var_data_rep = var_data_rep.replace("/", ".").replace("%","")
            if var.startswith("INFO/"):
                query_string = "%" + var
            else:
                query_string = "%INFO/" + var
            return {"query": query_string, "df": var_data_rep , "type": var_info["type"], "number": var_number}
        else:
            var_info = self.format_set[var.replace("FORMAT/","")]
            var_number = int(var_info["number"])
            if var_number > 1:
                var_data_rep =  ','.join([var + "." + str(v) for v in range(0,int(var_number))])
            else:
                var_data_rep = var
            var_data_rep = var_data_rep.replace("/", ".").replace("%","")
            query_string = "[%" + var.replace("FORMAT/","") + "]"
            return {"query": query_string, "df": var_data_rep , "type": var_info["type"], "number": var_number}
  
  def format_data_file_name(self,filename, x,y = None):
    filename = replace_all(filename,["bcf", "vcf","gz"], "").strip(".")
    if y == None:
        y = ""
    x, y  = replace_all(x, ["%","[","]"], ""), replace_all(y, ["%","[","]"], "")
    x, y  = x.replace("/", "-"), y.replace("/", "-")
    if y != "":
        y = "." + y
    return "{filename}.{x}{y}".format(**locals())

  def query(self, x, y=None, region=None, include=None):
    # Create analysis directory
    analysis_dir = replace_all(self.filename,["bcf", "vcf","gz"], "").strip(".")
    make_dir(analysis_dir)

    x = self.format_var_for_query(x)
    y = self.format_var_for_query(y)

    # Get dataframe rep to create header.
    if y["df"] == None:
        header = x["df"] + "\n"
    else:
        header = x["df"] + "," + y["df"] + "\n"

    if y["query"] == None:
        variables = x["query"]
    else:
        variables = x["query"] + "," + y["query"]

    q = shlex.split("bcftools query -f \"{variables}\\n\" {filename}".format(variables=variables, filename=self.filename))
    filename_pre = analysis_dir + "/" + self.format_data_file_name(self.filename, x["query"],y["query"])
    remove_file(filename_pre)
    with open(filename_pre + ".txt","w+") as out:
        out.write(header)
        Popen(q, stdout=out)

    if y["query"] == None:
        return filename_pre, x
    else:
        return filename_pre, x, y

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
