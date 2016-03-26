#! /usr/bin/env python
"""
usage:
  tb vcf2sql (bigquery|postgres|mysql|sqlite|datastore|dynamodb) [options] <vcf>

options:
  -h --help                   Show this screen.
  --version                   Show version.
  --filename=<filename>       SQLite Filename

  --db=<db>                   Database Name for MySQL, Postgres
  --user=<user>               User for MySQL, Postgres
  --password=<pw>             Password for MySQL, Postgres  
  --host=<host>               Host for MySQL, Postgres

  --table-name                Append a prefix to table names
  --vcf-version               Create a column indicating VCF version
  --simple                    Use Reduced Field Set
  --ann                    Parse snpeff fields

"""
from docopt import docopt
from utils.vcf import *
from subprocess import Popen, PIPE
import sys
import os
import re
import itertools
import gzip
from peewee import *
from pprint import pprint as pp
from signal import signal, SIGPIPE, SIG_DFL
from slugify import slugify
signal(SIGPIPE, SIG_DFL)


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

load_bigquery = r"""
#!/usr/bin/bash
# Load Bigquery

# [ ] Add instructions on getting set up...

gsutil cp {tsv_out} gs://andersen/annotation/{tsv_out}

# Load on bigquery
bq load --field_delimiter "\t" andersen-lab:Variation.{vcf_safe} gs://andersen/{tsv_out} {bigquery_schema}

"""


def bcftools_query(var_collection):
    query = []
    format_var_start = False
    for var, var_num, var_type, decription, group in var_collection:
        if group == "STANDARD":
            var = "%" + var
        elif group == "INFO":
            var = "%INFO/" + var
        elif group == "FORMAT":
            if format_var_start is False:
                var = "-->[%" + var
                format_var_start = True
            else:
                var = "%" + var
        query.append(var)
    query_in = r"\t".join(query) + r"\t%SAMPLE<--]\n"
    comm = Popen(["bcftools", "query", "-f", query_in,
                  v.filename], stdout=PIPE, stderr=PIPE)
    for line in comm.stdout:
        line = line.strip().split("-->")
        for sample_line in line[1].split("<--"):
            set_line = line[0] + sample_line
            set_line = zip(var_collection, set_line.split("\t"))
            out_line = []
            if sample_line:
                for q, o in set_line:
                    if q[2] == "Flag" and o == ".":
                        out_line.append("0")
                    elif q[2] == "Flag" and o != ".":
                        out_line.append("1")
                    elif o == ".":
                        out_line.append("")
                    else:
                        out_line.append(o)
                yield "\t".join(out_line)


def bigquery_schema(var_collection):
    var_set = []
    for var, var_num, var_type, decription, group in var_collection:
        try:
            var_num = int(var_num)
        except:
            var_num = 1
            var_type = "STRING"
        if var_type == "Flag":
            var_type = "BOOLEAN"
        if var_num > 1:
            var_type = "STRING"
        if group == "FORMAT":
            var = "F_" + var
        var_set.append([var, var_type.upper()])
    return var_set


def alert_script(script_name, script):
    """
        Prints script to user, notifies of location
    """
    with open(script_name, "w") as f:
        f.write(script)
    with indent(4):
        puts(colored.blue('\nA script named ' + script_name +
                          ' has been created. Run it to load the VCF when the the tsv is output.\n'))
    puts(colored.blue("<----- Script Contents ----->" +
                      script + "<--- End Script Contents --->"))

# Standard Columns
standard = (("CHROM", 1, "String", "Chromosome/Contig", "STANDARD"),
            ("POS", 1, "Integer", "Chromosome/Contig Position", "STANDARD"),
            ("_ID", 1, "String", "Variant ID", "STANDARD"),
            ("REF", 1, "String", "Reference Allele", "STANDARD"),
            ("ALT", 1, "String", "Alternative Alleles (list)", "STANDARD"),
            ("QUAL", 1, "Float", "Variant Quality", "STANDARD"),
            ("FILTER", 1, "String", "Variant Filters", "STANDARD"))


# Info
r_info = re.compile(r'''\#\#INFO=<
  ID=(?P<id>[^,]+),
  Number=(?P<number>-?\d+|\.|[AG]),
  Type=(?P<type>Integer|Float|Flag|Character|String),
  Description="(?P<desc>[^"]*)".*
  >''', re.VERBOSE)

# Format
r_format = re.compile(r'''\#\#FORMAT=<
  ID=(?P<id>.+),
  Number=(?P<number>-?\d+|\.|[AG]),
  Type=(?P<type>.+),
  Description="(?P<desc>.*)".*
  >''', re.VERBOSE)


try:
    os.remove("test.db")
except:
    pass

index_fields = ["CHROM", "POS", "SAMPLE", "FT", "GT", "_ID"]
simple_fields = ["FT", "GT", "TGT", "INDEL", "ANN"]
format_field_conversion = {"Integer": int, "String": str, "Float": float}


def filter_format_null(val):
    if val == ".":
        return None
    else:
        return val


def remove_missing_fields(val):
    ret_fields = []
    for k, v in val:
        if k in [x[0] for x in format_cols]:
            ret_fields.append((k, v))
    return ret_fields

def format_tgt(val, gt_dict):
    if val["GT"] == "./.":
        return "./."
    tgt = '/'.join([gt_dict[int(x)] for x in re.split("[\|/]", val["GT"])])
    return tgt

debug = None
if len(sys.argv) == 1:
    debug = ['vcf2sql', "test.vcf.gz"]

if __name__ == '__main__':
    args = docopt(__doc__,
                  argv=debug,
                  options_first=False)
    print args

    module_path = os.path.split(os.path.realpath(__file__))[0]
    v = vcf(args["<vcf>"])
    vcf_safe = v.filename.replace(".", "_")
    tsv_out = v.filename.replace("vcf", "tsv").replace(
        "bcf", "tsv").replace(".gz", "") + ".gz"

    info_cols = [map(autoconvert, list(x)) + ["INFO"]
                 for x in r_info.findall(v.raw_header)]
    format_cols = [map(autoconvert, list(x)) + ["FORMAT"]
                   for x in r_format.findall(v.raw_header)]

    if args["--simple"]:
        info_cols = [x for x in info_cols if x[0] in simple_fields]
        format_cols = [x for x in format_cols if x[0] in simple_fields]

    if args["bigquery"]:
        # Combine variable sets
        var_set = list(itertools.chain(standard, info_cols, format_cols)) + \
            [("SAMPLE", 1, "STRING", "Sample Name", "STANDARD")]
        # Generate bigquery schema
        bigquery_schema = ','.join([':'.join(x)
                                    for x in bigquery_schema(var_set)])

        script_filename = tsv_out.replace("tsv", "bigquery_load.sh")

        alert_script(script_filename, load_bigquery.format(**locals()))
        # Generate bcftools query
        with gzip.open(tsv_out, "wb") as f:
            for line in bcftools_query(var_set):
                f.write(line + "\n")
    elif args["sqlite"]:
        db = SqliteDatabase(args["--filename"])
    elif args["mysql"]:
        db = MySQLDatabase(args["--db"],
                           user = args["--user"],
                           password = args["--password"],
                           host = args["--host"])

    # Set table name
    tbl_name = slugify(args["<vcf>"])
    if args["--table-name"]:
        tbl_name = args["--table-name"]

    class vcf_table(Model):

        CHROM = CharField(index = True)
        POS = IntegerField(index = True)
        _ID = CharField(index = True, null = True)
        REF = CharField(index = True)
        ALT = CharField(index = True, null = True)
        QUAL = FloatField(null = True)
        FILTER = CharField(null = True)
        SAMPLE = CharField(index=True, help_text="Sample Name")
        TGT = CharField(index=True) 

        if args["--ann"]:
            allele = CharField(index=True)
            annotation = CharField(index=True)
            putative_impact = CharField(null=True)
            gene_name = CharField(index=True, null=True)
            gene_id = CharField(index=True, null=True)
            feature_type = CharField(null=True)
            feature_id = CharField(null=True)
            transcript_biotype = CharField(null=True)
            rank_total = CharField(null=True)
            hgvs_c = CharField(null=True)
            hgvs_p = CharField(null=True)
            cdna_position = CharField(null=True)
            cds_position = CharField(null=True)
            protein_position = CharField(null=True)
            distance_to_feature = CharField(null=True)
            errors = CharField(null=True)

        class Meta:
            database = db
            db_table = tbl_name

    #ann_field_types = [peewee_field_types[
    #    type(x)] for x in annotation._meta.sorted_fields[2:]]

    # Add in INFO Cols
    for i in info_cols:
        index_field = i[0] in index_fields
        if i[1] > 1:
            i[2] = "String"
        if i[2] == "String":
            CharField(index=index_field, null=True,
                      help_text=i[3]).add_to_class(vcf_table, i[0])
        elif i[2] == "Integer":
            IntegerField(index=index_field, null=True,
                         help_text=i[3]).add_to_class(vcf_table, i[0])
        elif i[2] == "Float":
            FloatField(index=index_field, null=True,
                       help_text=i[3]).add_to_class(vcf_table, i[0])

    for i in format_cols:
        if i[1] > 1:
            i[2] = "String"
        if i[2] == "String":
            index_field = i[0] in index_fields
            CharField(index=index_field, null=True,
                      help_text=i[3]).add_to_class(vcf_table, i[0])
        elif i[2] == "Integer":
            IntegerField(index=index_field, null=True,
                         help_text=i[3]).add_to_class(vcf_table, i[0])
        elif i[2] == "Float":
            FloatField(index=index_field, null=True,
                       help_text=i[3]).add_to_class(vcf_table, i[0])
    try:
        db.drop_tables([vcf_table], safe = True, cascade = True)
    except:
        pass
    db.create_tables([vcf_table])
    for loc in v:
        site_fields = {}
        annotation_record_set = []

        # Propogate standard fields
        for x in standard:
            if x[0] == "_ID":
                attr_name = "ID"
            else:
                attr_name = x[0]
            attr = getattr(loc, attr_name)
            if type(attr) == list:
                attr = ','.join(attr)
            site_fields[x[0]] = attr

        # Propogate INFO fields
        for x in info_cols:
            # Add Annotation Fields
            if x[0] == "ANN":
                annotation_record = {}
                for ann in loc.INFO["ANN"].split(","):
                    for field_name, field_type, field_value in zip(ann_fields, ann_field_types, ann.split("|")):
                        if not field_value:
                            annotation_record[field_name] = None
                        else:
                            annotation_record[
                                field_name] = field_type(field_value)
                    annotation_record_set.append(annotation_record)
            else:
                try:
                    attr = loc.INFO[x[0]]
                    if type(attr) == tuple or type(attr) == list:
                        attr = ','.join(map(str, attr))
                    site_fields[x[0]] = attr
                except:
                    site_fields[x[0]] = None
        with db.atomic():
            site_id = vcf_table(**site_fields)
            site_id.save()

            # Insert genotypes
            gt_dict = {0: loc.REF}
            alt_gt = dict(list(enumerate([0] + loc.ALT))[1:])
            gt_dict.update(alt_gt)
            loc_s = str(loc).strip().split("\t")
            format_str = loc_s[8].strip().split(":")
            # gt calls
            gt_set = [dict(remove_missing_fields(
                zip(format_str, map(filter_format_null, x.split(":"))))) for x in loc_s[9:]]
            # add sample names
            [x.update({"SAMPLE": sample, "TGT": format_tgt(x, gt_dict)})
             for sample, x in zip(v.samples, gt_set)]
            # Add site id to annotation records
            [x.update({"site": site_id}) for x in annotation_record_set]
            annotation.insert_many(annotation_record_set).execute()
            call.insert_many(gt_set).execute()
    
with indent(4):
    puts(colored.blue("\nDB Setup Complete\n"))
