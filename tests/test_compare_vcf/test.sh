#!/bin/bash
set -euo pipefail
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)"

export PATH=../../opt/jdk1.8.0_131/bin:$PATH
time ../../opt/miniconda2/bin/python ../../compare_vcf.py \
    --reference small.fa.gz --out_dir output \
    --vcfs call.vcf --true_vcf truth.vcf

TEST_FAIL=0
compare() {
    if [[ ($1 == *.gz) && ($2 == *.gz) ]];then
	if cmp -s <(zcat $1 | grep -v '##reference') <(zcat $2 | grep -v '##reference' );then
	    true
	else
	    TEST_FAIL=1
	fi
    elif cmp -s $1 $2; then
        true
    else    
        TEST_FAIL=1
    fi    

    if [[ TEST_FAIL -eq "1" ]]; then
        echo $1 $2 differs
    else    
        echo $1 $2 matches
    fi    
}
for i in expected/augmented*.vcf.gz;do
    compare $i output/$(basename $i)
done
if [[ $TEST_FAIL -ne 0 ]] ;then
    echo test fail
else    
    echo test pass
fi
