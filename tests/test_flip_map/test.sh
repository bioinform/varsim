#!/bin/bash
set -euo pipefail
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)"

time ../../flip_map.py flipped input.map 

TEST_FAIL=0
compare() {
    if [[ ($1 == *.gz) && ($2 == *.gz) ]];then
	if cmp -s <(zcat $1) <(zcat $2);then
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
for i in expected/*.map;do
    compare $i $(basename $i)
done
if [[ $TEST_FAIL -ne 0 ]] ;then
    echo test fail
else    
    echo test pass
fi
