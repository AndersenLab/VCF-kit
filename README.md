vcf-toolbox
===========

A package for graphing, QC, and comparing variant data from VCF files. Makes heavy use of bcftools.

### To Do

[ ] List available variables (INFO/General)

#### Common Options

[ ] __--title__ - Set plot options.
[ ] __--include__ - Use a custom filter.
[ ] __--region__ - Restrict analysis to a particular region
[ ] __--samples__ - Restrict analysis to particular sample(s)
[ ] __--facet__ - Facet by filter or variable.

#### Single Variable

[ ] plot <x> (float/int (histogram), categorical [bar])
[ ] --title (Set plot title)
[ ] --include= (allow use to use a custom filter).
[ ] --keep-data (save data)

#### Multiple Variables

[ ] plot <x> <y> (scatter plot | bar plot[freq x category])

#### Genome Distribution

[ ] plot genome <y> ; Plots a variable across the genome.
[ ] --bin-size (for binning); add warnings for large operations.

#### Quality Control

[ ] TSTV Ratio
[ ] --facet (by...variable)

#### Compare

[ ] plot compare
[ ] --qc-tstv; compare tstv ratios


#### Special

[ ] Generate summary report; incorporating a variety of different types of reports.
	- [ ] Allele frequency spectrum
	- [ ] Base Changes (Frequency)
	- [ ] Insertion/Deletion Lengths
	
