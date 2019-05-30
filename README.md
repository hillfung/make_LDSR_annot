# make_LDSR_annot (v: 2.4.6)
Collection of R functions to create annotations for analysis using stratified LDSC

## Getting started
Clone this repository and:
1. run `script_get_reference_file.sh`
2. source `make.LDSR.annot.R` in R

## Creating the annotation
To create an annotation, run the `make.LDSR.annot()`-function

This function requires at least 5 arguments:
- input=vector or data.frame (see below)
- out=output directory
- annot.name=name for the new annot
- template.dir=directory containing existing annotation files
- input.type=c("rs","id","name","bp","pos")

Other accepted arguments:
- sep=separator for the chromosome-basepair positions
- add.windows=add window-annotations (see below)
- GeneCode=reference file (see below)
- ID.extension=does the gene IDs contain extensions?
- gene.shoulder=add shoulder around genes (see below)

### input
The function can handle 5 types of input:
- a list of RSIDs (input.type='rs')
- a list of gene IDs (input.type='id')
- a list of gene names (input.type='name')
- a list of chromosome:basepair ranges (input.type='bp')
- a list of chromosome:basepair positions (input.type='pos')

For a list of CHR:BP ranges, the input must be a three-column data.frame with the following order:
- CHR
- START
- END
In all other cases the input should be either (1) a vector or (2) a single-column data.frame or matrix

**NOTE:** for a list of gene IDs/names, an additional file is required that indicates which basepair positions are spanned by the genes

This file can be created by running `script_get_reference_file.sh`

### add.windows
As indicated by Finucane et al (2015; 2018), it is sometimes proper to add so-called window-annotations around your main annotations. The argument accepts:
- an integer
- a vector of integers
- T(RUE)=add 100 and 500bp windows
- F(ALSE)

### GeneCode
This argument is used to indicate which **OBJECT** is the reference file that indicates the span of each gene. This object should contain the following 4 columns with these names:
- GENE_ID/GENE_NAME
- CHROMOSOME
- START
- END

### gene.shoulder
When using genes as input, it is often prudent to extend the genes with a certain range on **BOTH ENDS**. gene.shoulder accepts a single integer as input

## Output
The function produces 22 files called <annot.name>.[1-22].annot.gz
