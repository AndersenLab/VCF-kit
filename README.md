vcf-toolbox
===========

A package for graphing, QC, and comparing variant data from VCF files. Makes heavy use of bcftools.

### Planned Options


	Usage:
	  tb.py listvars <vcf>          
	  tb.py plot <vcf> <x> [<y>]      [options]
	  tb.py QC-Report <vcf>           [options]
	  tb.py compare <vcf>             [options]
	  tb.py -h | --help
	  tb.py --version

	Options:
	  -h --help                   Show this screen.
	  --version                   Show version.
	  --title=<title>             Set Custom plot titles.
	  --region=<region>           Restrict analysis to a particular region.
	  --include=<filter-expr>     Use a custom filtering string with bcftools.
	  --facet=<facet-var>         Facet analysis on a categorical variable.
	  --split-format              When plotting genotype FORMAT fields, facet by sample.



## Version 1.0 To Do List

#### General

- [X] List available variables (INFO/General)
- [X] Create/manage files within a created directory [create folder for vcf]
- [ ] Create single script for generating plots; append to it date/time, plot title, etc.

#### Common Options

- [ ] __--title__ - Set plot options.
- [ ] __--logx__ - Put x variable on log scale.
- [ ] __--logy__ - Put y variable on log scale.
- [ ] __--knitr__ - Create a knitr report; generate option.
- [ ] __--include__ - Use a custom filter.
- [ ] __--exclude__ - Filter Out Variants.
- [ ] __--region__ - Restrict analysis to a particular region
- [ ] __--samples__ - Restrict analysis to particular sample(s)
- [ ] __--facet__ - Facet by filter or categorical variable.

#### Listing Variables

	tb.py listvars <vcf>

- [X] List Variables, Types, etc.

#### Plotting Variables

	tb.py plot <vcf> <x> [<y>]

- [ ] plot individual variable
	- [X] Numeric
	- [ ] Categorical
    - [X] plot variant density across genome
		- [X] Genetic Labels (Mb, Kb, etc; Commas)
		- [X] bin-width
	- [ ] Multi-value fields
		- [ ] Reshaping
- [ ] Plot multiple variables
	- [ ] Numeric x Categorical
	- [ ] Numeric x Numeric (Scatterplot)
		- [ ] Labels (e.g. color by filter, another variable, tstv ratio)
	- [ ] Plot along chromosomal coordinates (e.g. DP by CHROM:POS)


#### Compare Variant Sets (e.g. diff callers)

	tb.py compare <x> [<y>] <vcf>

- [ ] Single Variable x Rate of Discordance
- [ ] POS --> Genomic Location and Rate of discordance
- [ ] Compare Two Variables
	- [ ] Numeric x Numeric
	- [ ] Numeric x Categorical
	- [ ] Markers for Concordance/Discordance
	- [ ] Facetting (esp. by filters)

#### Quality Control

	tb.py QC <vcflist>...

Outputs a QC report of one or more vcfs. If two vcfs are given, will compare them (and take considerably longer to run).

	- [ ] ts/tv ratio by chromosome:position; summary table
	- [ ] ts/tv ratio by depth
	- [ ] ts/tv ratio by quality
	- [ ] bcftools plots

#### Special

Odds and Ends

- [ ] Allele frequency spectrum
- [ ] Base Changes (Frequency)
- [ ] Insertion/Deletion Lengths
	
#### Future

* Plot on z axis?
* Multiple vcf comparison (Merge first?)
