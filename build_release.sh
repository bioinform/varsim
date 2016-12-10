#!/bin/bash

mvn package

version=$(git describe | sed 's/^v//')

tar -cv quickstart.sh VarSim.jar generate_small_test_ref.py varsim.py varsim_somatic.py | gzip > varsim-$version.tar.gz
