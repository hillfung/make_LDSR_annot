#!/bin/bash

## file obtained on 30-MAY-2019
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

## remove descriptives
zcat gencode.v19.annotation.gtf.gz | tail -n+6 > gencode.v19.annotation.gtf

## only keep "genes"
awk -F'\t' 'BEGIN{OFS="\t"};$3=="gene"{print $1,$4,$5,$7,$9}' < gencode.v19.annotation.gtf > gencode.v19.annotation.first_selected_columns.gtf
awk -F'; ' 'BEGIN{OFS="\t"}{print $1,$5}' < gencode.v19.annotation.first_selected_columns.gtf > gencode.v19.annotation.second_selected_columns.gtf
awk 'BEGIN{OFS="\t"};{print $1,$2,$3,$4,$6,$8}' < gencode.v19.annotation.second_selected_columns.gtf > gencode.v19.annotation.third_selected_columns.gtf

## remove "chr" from chromosome column
sed 's!^chr!!g' < gencode.v19.annotation.third_selected_columns.gtf > gencode.v19.annotation.third_selected_columns.numeric_chr.gtf

## remove quotes around gene IDs and gene names
sed 's!"!!g' < gencode.v19.annotation.third_selected_columns.numeric_chr.gtf > gencode.v19.annotation.third_selected_columns.numeric_chr.no_quotes.gtf

## add header
echo -e "CHR\tSTART\tEND\tSTRAND\tGENE_ID\tGENE_NAME" | cat - gencode.v19.annotation.third_selected_columns.numeric_chr.no_quotes.gtf > gencode_v19_genes_only.tab
