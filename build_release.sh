#!/bin/bash

version=$(git describe | sed 's/^v//')

mvn versions:set -DgenerateBackupPoms=false -DnewVersion=$version

mvn package

tar -cv quickstart.sh VarSim.jar generate_small_test_ref.py varsim.py varsim_somatic.py | gzip > varsim-$version.tar.gz
