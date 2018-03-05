#!/bin/bash

cat $1 | grep "^#"; cat $@ | grep -v "^#" | sort -k1,1V -k2,2n -k4,4d -k5,5d
