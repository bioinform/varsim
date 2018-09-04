#!/usr/bin/env bash

set -ex

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OPT_DIR=${DIR}/opt

mkdir -p ${OPT_DIR}
pushd ${OPT_DIR}

wget --no-check-certificate --no-cookies --header "Cookie: oraclelicense=accept-securebackup-cookie" http://download.oracle.com/otn-pub/java/jdk/8u131-b11/d54c1d3a095b4ff2b6607d096fa80163/jdk-8u131-linux-x64.tar.gz
tar -zxvf jdk-8u131-linux-x64.tar.gz
rm jdk-8u131-linux-x64.tar.gz

export PATH=${OPT_DIR}/jdk1.8.0_131/bin:${PATH}

PYTHON_DIR=${OPT_DIR}/miniconda2 && rm -rf ${PYTHON_DIR}
CONDA=Miniconda2-latest-Linux-x86_64.sh
wget -q https://repo.continuum.io/miniconda/${CONDA}\
    && sh ${CONDA} -b -p ${PYTHON_DIR}\
    && ${PYTHON_DIR}/bin/pip install -I pysam==0.15.0\
    && ${PYTHON_DIR}/bin/pip install -I pyvcf==0.6.8\
    && ${PYTHON_DIR}/bin/conda install --yes -c daler pybedtools=0.7.2 bedtools=2.25.0 \
    && ${PYTHON_DIR}/bin/pip install -I pandas==0.23.3\
    && ${PYTHON_DIR}/bin/pip install -I numpy==1.15.0\
    && ${PYTHON_DIR}/bin/pip install -I scipy==1.1.0\
    && rm -f ${CONDA}

wget -O- http://mirrors.sonic.net/apache/maven/maven-3/3.5.4/binaries/apache-maven-3.5.4-bin.tar.gz | tar zxvf -
wget -O- https://www.apache.org/dist/ant/binaries/apache-ant-1.9.13-bin.tar.gz | tar zxvf -
popd

export JAVA_HOME=${OPT_DIR}/jdk1.8.0_131
version=$(git describe | sed 's/^v//')
${OPT_DIR}/apache-maven-3.5.4/bin/mvn versions:set -DgenerateBackupPoms=false -DnewVersion=$version
${OPT_DIR}/apache-maven-3.5.4/bin/mvn package

git submodule init
git submodule update
pushd rtg-tools
rm -rf rtg-tools-*
${OPT_DIR}/apache-ant-1.9.13/bin/ant zip-nojre
unzip dist/rtg-tools*.zip
cp rtg-tools*/RTG.jar $DIR
popd
