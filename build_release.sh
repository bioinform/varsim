#!/bin/bash
set -e

version=$(git describe | sed 's/^v//')
cwd=`pwd`

git checkout -- pom.xml
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

git checkout -- pom.xml
tar -cv pom.xml tests/quickstart_test/quickstart.sh build.sh utils.py compare_vcf.py RTG.jar VarSim.jar __init__.py liftover_restricted_vcf_map.py generate_small_test_ref.py varsim.py varsim_multi.py varsim_somatic.py varsim_validator.py | gzip > varsim-$version.tar.gz
