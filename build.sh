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
    && ${PYTHON_DIR}/bin/pip install -I pysam\
    && ${PYTHON_DIR}/bin/pip install -I pyvcf\
    && ${PYTHON_DIR}/bin/conda install --yes -c daler pybedtools bedtools \
    && ${PYTHON_DIR}/bin/pip install -I pandas\
    && ${PYTHON_DIR}/bin/pip install -I numpy\
    && ${PYTHON_DIR}/bin/pip install -I scipy\
    && rm -f ${CONDA}
