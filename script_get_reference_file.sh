#!/bin/bash

printf '*%.0s' {1..96}
echo -e '
*
* Script to create reference file to be used in making LDSC annotations (v1.1.0)
* (C) 2019 Hill F. Ip
* Department of Biological Psychology, Vrije Universiteit Amsterdam / Netherlands Twin Register
* GNU General Public License v3
*'
printf '*%.0s' {1..96}
echo -e '\n'

BEGINTIME=$(date "+%s.%N")

echo -e "Script started at $(date --date @$(echo ${BEGINTIME} | cut -f1 -d'.') "+%a %d-%b-%Y %T %Z")\n"

## file obtained on 31-MAY-2019
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

## remove descriptives
echo -en "Unzipping gzipped file..."
zcat gencode.v19.annotation.gtf.gz | tail -n+6 > gencode.v19.annotation.gtf
echo "done"

## only keep "genes" + drop unused columns
echo -en "Removing anything that are not genes..."
awk -F'\t' 'BEGIN{OFS="\t"};$3=="gene"{print $1,$4,$5,$7,$9}' < gencode.v19.annotation.gtf > gencode.v19.annotation.genes_only.gtf
echo "done"

## extract gene IDs and gene names from descriptives
echo -en "Extracting gene IDs and gene names..."
awk -F'; ' 'BEGIN{OFS="\t"}{print $1,$3,$4,$5}' < gencode.v19.annotation.genes_only.gtf > gencode.v19.annotation.genes_only.geneIDs_geneNames.gtf
echo "done"

## drop unused columns (i.e. columns of "gene_id", "gene_type", "gene_status", and "gene_name")
echo -en "Dropping unused columns"
awk 'BEGIN{OFS="\t"};{print $1,$2,$3,$4,$6,$8,$10,$12}' < gencode.v19.annotation.genes_only.geneIDs_geneNames.gtf > gencode.v19.annotation.genes_only.geneIDs_geneNames.selected_columns.gtf
echo "done"

## remove "chr" from chromosome column
echo -en "Removing \"chr\" from the chromosomee column..."
sed 's!^chr!!g' < gencode.v19.annotation.genes_only.geneIDs_geneNames.selected_columns.gtf > gencode.v19.annotation.genes_only.geneIDs_geneNames.selected_columns.numeric_chr.gtf
echo "done"

## remove quotes
echo -en "Removing quotes from gene IDs and gene names..."
sed 's!"!!g' < gencode.v19.annotation.genes_only.geneIDs_geneNames.selected_columns.numeric_chr.gtf > gencode.v19.annotation.genes_only.geneIDs_geneNames.selected_columns.numeric_chr.no_quotes.gtf
echo "done"

## remove duplicates based on gene names (gene IDs are all unique)
echo -e "Removing duplicate entries based on gene names\nNOTE: gene IDs are all unique at this point\n"
cut -f8 gencode.v19.annotation.genes_only.geneIDs_geneNames.selected_columns.numeric_chr.no_quotes.gtf | sort | uniq -c | awk '$1!=1{print $2}' > duplicated_geneNames.list
awk 'BEGIN{OFS="\t"}
FNR==NR{a[$1]=$1;next}
FNR!=NR{if($8 in a){print $0 > "gencode.v19.annotation.genes_only.geneIDs_geneNames.selected_columns.numeric_chr.no_quotes.duplicated_geneNames.gtf"}else{print $1,$2,$3,$4,$5,$8 > "gencode_v19_genes_only.no_header.tab"}}
' duplicated_geneNames.list gencode.v19.annotation.genes_only.geneIDs_geneNames.selected_columns.numeric_chr.no_quotes.gtf
while read GENE_NAME;do
awk -v GENE_NAME=${GENE_NAME} '$8==GENE_NAME' < gencode.v19.annotation.genes_only.geneIDs_geneNames.selected_columns.numeric_chr.no_quotes.duplicated_geneNames.gtf > dum
if [[ $(cut -f1 dum | sort | uniq | wc -l) != 1 ]];then
	echo -e "  WARNING: \"${GENE_NAME}\" was found on different chromosomes and excluded from the file"
	echo "${GENE_NAME}" >> excluded.list
else
	MAX_BEGIN=$(cut -f2 dum | sort -g | tail -n1)
	MIN_END=$(cut -f3 dum | sort -g | head -n1)
	if [[ ${MAX_BEGIN} -gt ${MIN_END} ]];then
		echo -e "  WARNING: \"${GENE_NAME}\" contains non-overlapping entries and was excluded from the file"
		echo "${GENE_NAME}" >> excluded.list
	else
		awk 'BEGIN{OFS="\t"};{print $1,$2,$3,$4,$5,$8}' < dum >> gencode_v19_genes_only.no_header.tab
	fi
fi
rm dum
done < duplicated_geneNames.list

## add header
echo -en "\nAdding header to output file"
echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tGENE_ID\tGENE_NAME" | cat - gencode_v19_genes_only.no_header.tab > gencode_v19_genes_only.tab
echo "done"

if [[ $(find . -maxdepth 1 -name "gencode_v19_genes_only.tab" -type f 2> /dev/null | wc -l) == 1 ]];then
	echo -e "Succesfully created reference file for $(( $(wc -l < gencode_v19_genes_only.tab) - 1 )) entries\nReference file is saved as \"gencode_v19_genes_only.tab\""

	ENDTIME=$(date "+%s.%N")
	TOTAL_TIME=$(date -u -d "0 $ENDTIME sec - ${BEGINTIME} sec" "+%H:%M:%S.%3N")
	echo -e "\nScript finished at $(date --date @$(echo ${ENDTIME} | cut -f1 -d'.') "+%a %d-%b-%Y %T %Z")\nElapsed time: ${TOTAL_TIME}"
else
	echo "ERROR: something happened and the script did not properly make the reference file"
	exit 1
fi

