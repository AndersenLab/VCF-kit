#! /usr/bin/env python
"""
usage:
  vk vcf2sql (postgres|mysql|sqlite) [options] <vcf>

options:
  -h --help                   Show this screen.
  --version                   Show version.

  --db=<db>                   Database Name for MySQL, Postgres OR filename for sqlite
  --user=<user>               User for MySQL, Postgres
  --password=<pw>             Password for MySQL, Postgres  
  --host=<host>               Host for MySQL, Postgres
  --compress                  Compress GT data.

  --table-name=<table-name>   Append a prefix to table names
  --vcf-version               Create a column indicating VCF version
  --simple                    Use Reduced Field Set
  --ANN                       Parse snpeff fields
  --modifier                  Include modifier annotation records.
  --print                     Print CSV to stdout instead of loading.

"""
from docopt import docopt
from utils.vcf import * 
import sys
import os
import re
from peewee import *
from signal import signal, SIGPIPE, SIG_DFL
from slugify import slugify
import datetime
import copy
import tempfile
import csv
import json
import zlib, cPickle, base64
from playhouse.csv_loader import load_csv
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

index_fields = ["CHROM", "POS", "SAMPLE", "FT", "_ID", "REF", "ALT", "gene_name", "gene_id", "putative_impact", "feature_type"]
simple_fields = ["FT", "GT", "TGT", "INDEL", "ANN"]
format_field_conversion = {"Integer": int, "String": str, "Float": float}

ann_fields = ["allele", 
              "annotation",
              "putative_impact",
              "gene_name",
              "gene_id",
              "feature_type",
              "feature_id",
              "transcript_biotype",
              "rank_total",
              "hgvs_c",
              "hgvs_p",
              "cdna_position",
              "cds_position",
              "protein_position",
              "distance_to_feature",
              "errors"]


ann_field_types = [str]*14


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
        return None
    tgt = '/'.join([gt_dict[int(x)] for x in re.split("[\|/]", val["GT"])])
    return tgt

debug = None
if len(sys.argv) == 1:
    debug = ['vcf2sql', "test.vcf.gz"]

if __name__ == '__main__':
    args = docopt(__doc__,
                  argv=debug,
                  options_first=False)

    timestamp = datetime.datetime.now()

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

    if args["sqlite"]:
        db = SqliteDatabase(args["--db"])
    elif args["mysql"]:
        db = MySQLDatabase(args["--db"],
                           user = args["--user"],
                           password = args["--password"],
                           host = args["--host"])

    if not args["--print"]:
        db.connect()
        # Set table name
        tbl_name = slugify(args["<vcf>"])
        if args["--table-name"]:
            tbl_name = args["--table-name"]


    class vcf_table(Model):
        Variant = IntegerField(index = True)
        CHROM = CharField(index = True)
        POS = IntegerField(index = True)
        _ID = CharField(index = True, null = True)
        REF = CharField(index = True)
        ALT = CharField(index = True, null = True)
        QUAL = FloatField(null = True)
        FILTER = CharField(null = True)
        GT = TextField(null = True)

        if args["--ANN"]:
            allele = CharField(index=True, null = True)
            annotation = CharField(index=True, null = True)
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
        if not args["--print"]:
            class Meta:
                database = db
                db_table = tbl_name

    # Add in INFO Cols
    for i in info_cols:
        index_field = i[0] in index_fields
        if i[0] == "ANN":
            break
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
        elif i[2] == "Flag":
            BooleanField(null = True).add_to_class(vcf_table, i[0])

    # Output Schema
    with open("schema.sql", "w") as f:
        f.write("".join(vcf_table.sqlall()))

    db.create_tables([vcf_table], safe = True)
    c = 0
    loc_set = [] # 
    insert_set = [] # 
    csv_writer = csv.writer(sys.stdout, quoting=csv.QUOTE_MINIMAL)
    for variant_index, loc in enumerate(v):
        c += 1
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

        annotation_record_set = []
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
                            annotation_record[field_name] = field_type(field_value)
                    if ((annotation_record["putative_impact"] != "MODIFIER") or args["--modifier"]):
                        annotation_record_set.append(copy.copy(annotation_record))
            else:
                try:
                    attr = loc.INFO[x[0]]
                    if type(attr) == tuple or type(attr) == list:
                        attr = ','.join(map(str, attr))
                    site_fields[x[0]] = attr
                except:
                    site_fields[x[0]] = None

        # Insert genotypes
        gt_dict = {0: loc.REF}
        alt_gt = dict(list(enumerate([0] + loc.ALT))[1:])
        gt_dict.update(alt_gt)
        loc_s = str(loc).strip().split("\t")
        format_str = loc_s[8].strip().split(":")

        # Process genotype calls
        # gt calls
        gt_set = [dict(remove_missing_fields(
            zip(format_str, map(filter_format_null, x.split(":"))))) for x in loc_s[9:]]

        # add sample names
        [x.update({"SAMPLE": sample, "TGT": format_tgt(x, gt_dict)})
         for sample, x in zip(v.samples, gt_set)]

        # Insert Data
        if args["--ANN"]:
            record_inserted = False # Ensure at least one record inserted
            for i in annotation_record_set:
                loc_set.append(copy.copy(i))
                record_inserted = True
            if not record_inserted:
                loc_set.append({})
        else:
            loc_set.append({})
        
        field_names = vcf_table._meta.sorted_field_names[1:]
        for rec in loc_set:
            rec["Variant"] = variant_index
            rec["CHROM"] = loc.CHROM
            rec["POS"] = loc.POS
            rec["_ID"] = loc.ID
            rec["REF"] = loc.REF
            rec["ALT"] = '|'.join(loc.ALT)
            rec["QUAL"] = loc.QUAL
            rec["FILTER"] = loc.FILTER
            if args["--compress"]:
                rec["GT"] = base64.b64encode(zlib.compress(cPickle.dumps(gt_set).encode("utf-8")))
            else:
                rec["GT"] = json.dumps(gt_set)

            
            for k in field_names:
                if k not in rec:
                    rec[k] = ""
                elif rec[k] is None:
                    rec[k] = ""
            if args["--print"]:
                csv_writer.writerow(["0"] + [str(rec[x]) for x in field_names])
            else:
                insert_set.append(rec)
            loc_set = []
        if c % 15000 == 0 and not args["--print"]:
            with db.atomic():
                vcf_table.insert_many(insert_set).execute()
            print("Inserted {c} records".format(c=c))
            insert_set = []
    if not args["--print"]:
        with db.atomic():
                vcf_table.insert_many(insert_set).execute()

with indent(4):
    puts_err(colored.blue("\nDB Setup Complete\n"))
