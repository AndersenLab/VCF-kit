vcf-toolbox
===========

A package for graphing, QC, and comparing variant data from VCF files. Makes heavy use of bcftools.

### To Do

- [X] List available variables (INFO/General)
- [X] Create/manage files within a created directory [create folder for vcf]

#### Common Options

- [ ] __--title__ - Set plot options.
- [ ] __--include__ - Use a custom filter.
- [ ] __--exclude__ - Filter Out Variants.
- [ ] __--region__ - Restrict analysis to a particular region
- [ ] __--samples__ - Restrict analysis to particular sample(s)
- [ ] __--facet__ - Facet by filter or categorical variable.

#### Single Variable

- [ ] plot <x> (float/int (histogram), categorical [bar])
	- [X] Integer
	- [ ] Float
	- [ ] Categorical
    - [X] plot variant density across genome
		- [X] Genetic Labels (Mb, Kb, etc; Commas)
		- [X] bin-width
	- [ ] Multi-value fields

#### Multiple Variables

- [ ] plot <x> <y> (scatter plot | bar plot[freq x category])
	- [ ] - Plot variable across genome (avg; binned)

#### Quality Control

- [ ] TSTV Ratio
	- [ ] - Binned; Across genome; By sample?
- [ ] __--facet__ (by...variable)

#### Compare

- [ ] plot compare
- [ ] __--qc-tstv__ - compare tstv ratios


#### Special

- [ ] Generate summary report; incorporating a variety of different types of reports.
	- [ ] Allele frequency spectrum
	- [ ] Base Changes (Frequency)
	- [ ] Insertion/Deletion Lengths
	
