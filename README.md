# make_LDSR_annot (v: 2.4.6)
Collection of R functions to create annotations for analysis using stratified LDSC

There are two ways to create a new annotation:
1. specifying a set of SNPs to include
2. collapse existing annotations

## Requirements
The function requires the `data.table`-package to be installed in R

## Getting started
Clone this repository and:
1. run `script_get_reference_file.sh` (optional)
2. source `make.LDSR.annot.R` in R

## 1. Creating a new annotation by specifying a set of SNPs
To create a new annotation, run the `make.LDSR.annot()`-function

This function requires the following 5 arguments:
- input=vector or data.frame (see below)
- out=output directory
- annot.name=name for the new annotation
- template.dir=directory containing existing annotation files (see below)
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

**NOTE:** for a list of gene IDs/names, an additional file is required that indicates which basepair positions are spanned by the genes (see below). This file can be created by running `script_get_reference_file.sh`

### template.dir
Use this argument to indicate where the templates are located. For example, `template.dir="/home/user/my_annotations/foo"` results in the function using `foo.[1-22].annot.gz` as templates for the new annotation.

### add.windows
As indicated by [Finucane et al (2015)](https://doi.org/10.1038/ng.3404), it is sometimes proper to add additional so-called window-annotations next to your main annotation. The window annotation includes all SNPs that are part of the main annotation and every SNP that is within a specific BP-window from these SNPs. The argument accepts one of the following values:
- an integer
- a vector of integers
- T(RUE)=add 100 and 500bp windows (default)
- F(ALSE)=do not add windows

### GeneCode
This argument is used to indicate which **OBJECT** is the reference file that indicates the span of each gene. This object should contain the following 4 columns with these names:
- GENE_ID/GENE_NAME
- CHROMOSOME
- START
- END

**NOTE:** the object may contain more columns than required (e.g. both GENE_NAME and GENE_ID) as long as the column names of the relevant columns match exactly with what is shown above

### gene.shoulder
When using genes as input, it is often prudent to extend the genes with a certain range (see [Finucane et al, 2018](https://doi-org.vu-nl.idm.oclc.org/10.1038/s41588-018-0081-4)). `gene.shoulder` accepts a single integer as input.

**Note:** this is different from specifying windows: `add.windows=100` will return a file with two annotations, the main annotation and a 100-BP window; whereas `gene.shoulder=100` returns one annotation where the ranges as specified by `GeneCode` are extended with 100 BP **on both sides** of the gene. 

**NOTE:** if both `add.windows` and `gene.shoulder` are specified, gene ranges will first be extended and windows will be made based on the extended ranges

## 2. collapse existing annotations
Use the `collapse.annots()`-function to combine multiple annotations into a single annotation. You can choose to do one of two things:
- **intersection**=new annotation will only include SNPs that are **common to all annotations**
- **union**=new annotation will contain all SNPs that is **part of at least one annotation**

**NOTE:** the function has only been tested with files that **contain only one annotation**. The function may not work propberly with files that contain windows and/or consists of multiple annotations (e.g. the "baseline"-annotations).

The function requires 4 arguments:
- annot.dir=directory containing the annotations to collapse
- out=output directory
- annot.name=name for the new annotation
- collapse.type= "intersection" or "union" (case-sensitive)

Optional argument:
- add.windows (see above)

### annot.dir
The value supplied to this argument should be a (vector of) strings. For example: `annot.dir=c("/home/user/my_annotations/foo","/home/user/my_annotations/bar")` indicates that `foo.[1-22].annot.gz` and `bar.[1-22].annot.gz` will be included in the new annotation.

Optionally, the value may be a directory and all annotations inside the annotation will be collapsed into the new annotation.

## Output
Both functions produce 22 files called <annot.name>.[1-22].annot.gz
