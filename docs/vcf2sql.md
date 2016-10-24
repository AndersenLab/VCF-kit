# Overview


### vcf2sql

Prepares a VCF for import into an sql database. Bigquery, Postgres, Mysql, and sqlite will be supported. Data is loaded denormalized and
indices are added to enable easy querying.

```
  vk vcf2sql (bigquery|postgres|mysql|sqlite) <vcf>
```

* [X] Bigquery
* [X] Mysql
* [X] Postgres
* [X] Sqlite
* [X] Automatically load
* [X] Support for SNPeff annotations
* [ ] Reorder columns
* [ ] Support for multi-column types
* [ ] Break load job for bigquery into multiple files if filesize > 4GB
* [ ] Output TGT (bases), allele 1, allele 2 for genotypes.

### vcf2bigquery

* [ ] Split bigquery into its own tool
