#!/bin/bash

version=$(git describe | sed 's/^v//')

mvn versions:set -DgenerateBackupPoms=false -DnewVersion=$version

mvn package

tar -cv quickstart.sh VarSim.jar __init__.py liftover_restricted_vcf_map.py generate_small_test_ref.py varsim.py varsim_somatic.py | gzip > varsim-$version.tar.gz
