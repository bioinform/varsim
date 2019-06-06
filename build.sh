#!/usr/bin/env bash

set -ex

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OPT_DIR=${DIR}/opt

mkdir -p ${OPT_DIR}
pushd ${OPT_DIR}


JDK8_DIR="${OPT_DIR}/jdk1.8.0_131"
if [ ! -d ${JDK8_DIR} ];then
wget -O- --no-check-certificate --no-cookies --header "Cookie: oraclelicense=accept-securebackup-cookie" \
    http://download.oracle.com/otn-pub/java/jdk/8u131-b11/d54c1d3a095b4ff2b6607d096fa80163/jdk-8u131-linux-x64.tar.gz | tar -zxvf -
fi

export PATH=${JDK8_DIR}/bin:${PATH}

PYTHON_DIR=${OPT_DIR}/miniconda2
CONDA=Miniconda2-latest-Linux-x86_64.sh
if [[ ! -d ${PYTHON_DIR} ]]; then
wget -q https://repo.continuum.io/miniconda/${CONDA}\
    && sh ${CONDA} -b -p ${PYTHON_DIR}\
    && ${PYTHON_DIR}/bin/python ${PYTHON_DIR}/bin/pip install pysam==0.15.0\
    && ${PYTHON_DIR}/bin/python ${PYTHON_DIR}/bin/pip install pyvcf==0.6.8\
    && ${PYTHON_DIR}/bin/python ${PYTHON_DIR}/bin/conda install --yes -c bioconda pybedtools=0.8.0 bedtools=2.25.0 \
    && ${PYTHON_DIR}/bin/python ${PYTHON_DIR}/bin/pip install scipy==1.1.0\
    && rm -f ${CONDA}
fi

MAVEN_DIR=${OPT_DIR}/apache-maven-3.5.4
if [[ ! -d ${MAVEN_DIR} ]]; then
wget -O- http://mirrors.sonic.net/apache/maven/maven-3/3.5.4/binaries/apache-maven-3.5.4-bin.tar.gz | tar zxvf -
fi
ANT_DIR=${OPT_DIR}/apache-ant-1.9.13
if [[ ! -d ${ANT_DIR} ]]; then
wget -O- https://www.apache.org/dist/ant/binaries/apache-ant-1.9.13-bin.tar.gz | tar zxvf -
fi

BZIP_DIR=${OPT_DIR}/bzip2-1.0.6
if [[ ! -d ${BZIP_DIR} ]]; then
    wget -O- https://www.sourceware.org/pub/bzip2/bzip2-1.0.6.tar.gz | tar zxvf -
    pushd ${BZIP_DIR}
    make install PREFIX=${BZIP_DIR}_install CFLAGS=" -fPIC"
    popd
fi

# Download samtools and index reference
samtools_version="1.9"
SAMTOOLS_DIR=${OPT_DIR}/samtools-${samtools_version}
if [[ ! -d $SAMTOOLS_DIR ]]; then
    wget -O- https://github.com/samtools/samtools/releases/download/$samtools_version/samtools-$samtools_version.tar.bz2 | tar xfj -
    pushd ${SAMTOOLS_DIR}
    C_INCLUDE_PATH_CURRENT=$C_INCLUDE_PATH
    export C_INCLUDE_PATH=${BZIP_DIR}_install/include:${C_INCLUDE_PATH}
    make install prefix=${SAMTOOLS_DIR}_install LDFLAGS="-L${BZIP_DIR}_install/lib"
    pushd htslib-1.9
    make install prefix=../../htslib-1.9_install LDFLAGS="-L${BZIP_DIR}_install/lib"
    popd
    popd
    export C_INCLUDE_PATH=$C_INCLUDE_PATH_CURRENT
fi

# Download ART
ART_DIR=${OPT_DIR}/ART
if [[ ! -d ${ART_DIR} ]]; then
    mkdir -p ${ART_DIR}
    pushd ${ART_DIR}
    wget -O- http://www.niehs.nih.gov/research/resources/assets/docs/artbinvanillaicecream031114linux64tgz.tgz | tar xfz -
    popd
fi

popd
export JAVA_HOME=${OPT_DIR}/jdk1.8.0_131
version=$(git describe | sed 's/^v//')
${MAVEN_DIR}/bin/mvn versions:set -DgenerateBackupPoms=false -DnewVersion=$version
${MAVEN_DIR}/bin/mvn package

git submodule init
git submodule update
pushd rtg-tools
rm -rf rtg-tools-*
${ANT_DIR}/bin/ant zip-nojre
unzip dist/rtg-tools*.zip
cp rtg-tools*/RTG.jar $DIR
popd
