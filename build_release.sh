#!/bin/bash
set -e

version=$(git describe | sed 's/^v//')
cwd=`pwd`

mvn versions:set -DgenerateBackupPoms=false -DnewVersion=$version

mvn package

git submodule init
git submodule update
pushd rtg-tools
rm -rf rtg-tools-*
ant zip-nojre
unzip dist/rtg-tools*.zip
cp rtg-tools*/RTG.jar $cwd
popd

tar -cv RTG.jar quickstart.sh VarSim.jar __init__.py liftover_restricted_vcf_map.py generate_small_test_ref.py varsim.py varsim_multi.py varsim_somatic.py varsim_validator.py | gzip > varsim-$version.tar.gz
