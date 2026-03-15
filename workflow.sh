#!/bin/bash

INPUT=$1
OUTPUT=$2

zcat $INPUT | grep -v "^#"| grep -P '^(chr[0-9]+|chr[ZMYX])\t' > ${OUTPUT}no_headers.vcf

NOHEADERS=${OUTPUT}no_headers.vcf


cat $NOHEADERS | cut -f1-6 > ${OUTPUT}part1_core.tsv

cat $NOHEADERS | grep -o -E "DP=[0-9]+" | sed 's/DP=//' > ${OUTPUT}part2_dp.tsv

cat $NOHEADERS | awk '{if($8 ~ /INDEL/) print "INDEL"; else print "SNP"}' > ${OUTPUT}part3_type.tsv

wc -l ${OUTPUT}part1_core.tsv ${OUTPUT}part2_dp.tsv ${OUTPUT}part3_type.tsv

printf "CHROM\tPOS\tID\tREF\tALT\tQUAL\tDP\tTYPE\n" > "${OUTPUT}.tsv"

paste ${OUTPUT}part1_core.tsv ${OUTPUT}part2_dp.tsv ${OUTPUT}part3_type.tsv >> "${OUTPUT}.tsv"

echo "Output saved to ${OUTPUT}.tsv"
wc -l "${OUTPUT}.tsv"