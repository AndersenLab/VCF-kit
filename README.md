vcf-toolbox
===========

A package for graphing, QC, and comparing variant data from VCF files. Makes heavy use of bcftools.

### To Do

- [X] List available variables (INFO/General)
- [X] Create/manage files within a created directory [create folder for vcf]

#### Common Options

- [ ] __--title__ - Set plot options.
- [ ] __--include__ - Use a custom filter.
- [ ] __--region__ - Restrict analysis to a particular region
- [ ] __--samples__ - Restrict analysis to particular sample(s)
- [ ] __--facet__ - Facet by filter or categorical variable.

#### Single Variable

- [ ] plot <x> (float/int (histogram), categorical [bar])
	- [X] Integer
	- [ ] Float
	- [ ] Categorical
	- [X] (POS) plot
		- [ ] Genetic Labels (Mb, Kb, etc; Commas)
	- [ ] Multi-value fields

#### Multiple Variables

- [ ] plot <x> <y> (scatter plot | bar plot[freq x category])

#### Genome Distribution

- [ ] plot genome <y> ; Plots a variable across the genome.
- [ ] __--bin-size__ (for binning); add warnings for large operations.

#### Quality Control

- [ ] TSTV Ratio
- [ ] __--facet__ (by...variable)

#### Compare

- [ ] plot compare
- [ ] __--qc-tstv__ - compare tstv ratios


#### Special

- [ ] Generate summary report; incorporating a variety of different types of reports.
	- [ ] Allele frequency spectrum
	- [ ] Base Changes (Frequency)
	- [ ] Insertion/Deletion Lengths
	
