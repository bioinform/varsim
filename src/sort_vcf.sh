#!/bin/bash
set -e

for vcf in $@;do
    if [ ! -f $vcf ]; then
	exit 1
    fi
done    

if [[ $1 == *.gz ]]; then
    zcat $1 
else
    cat $1 
fi | grep "^#"

for vcf in $@;do
    if [[ $vcf == *.gz ]]; then
        zcat $vcf
    else
        cat $vcf
    fi
done | grep -v "^#" | sort -k1,1V -k2,2n -k4,4d -k5,5d 
