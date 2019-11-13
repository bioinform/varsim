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
CONDA=Miniconda2-4.6.14-Linux-x86_64.sh
samtools_version="1.9"
SAMTOOLS_DIR=${OPT_DIR}/samtools-${samtools_version}_install
HTSLIB_DIR=${OPT_DIR}/htslib-1.9_install

if [[ ! -d ${PYTHON_DIR} ]]; then
wget -q https://repo.continuum.io/miniconda/${CONDA}\
    && sh ${CONDA} -b -p ${PYTHON_DIR}\
    && ${PYTHON_DIR}/bin/python ${PYTHON_DIR}/bin/pip install pysam==0.15.0\
    && ${PYTHON_DIR}/bin/python ${PYTHON_DIR}/bin/pip install pyvcf==0.6.8\
    && ${PYTHON_DIR}/bin/python ${PYTHON_DIR}/bin/conda install --yes -c bioconda pybedtools=0.8.0 bedtools=2.25.0 \
    && ${PYTHON_DIR}/bin/python ${PYTHON_DIR}/bin/pip install scipy==1.1.0\
    && ${PYTHON_DIR}/bin/python ${PYTHON_DIR}/bin/conda install --yes -c bioconda samtools=${samtools_version} htslib=1.9\
    && rm -f ${CONDA}

    if [[ ! -d $SAMTOOLS_DIR ]]; then
        mkdir -p ${SAMTOOLS_DIR}/bin/
        for i in samtools bgzip tabix;do
            ln -sf ${PYTHON_DIR}/bin/${i} ${SAMTOOLS_DIR}/bin/
        done
    fi

    if [[ ! -d $HTSLIB_DIR ]]; then
        mkdir -p ${HTSLIB_DIR}/bin/
        for i in bgzip htsfile tabix;do
            ln -sf ${PYTHON_DIR}/bin/${i} ${HTSLIB_DIR}/bin/
        done
    fi
fi

MAVEN_DIR=${OPT_DIR}/apache-maven-3.5.4
if [[ ! -d ${MAVEN_DIR} ]]; then
wget -O- http://mirrors.sonic.net/apache/maven/maven-3/3.5.4/binaries/apache-maven-3.5.4-bin.tar.gz | tar zxvf -
fi
ANT_DIR=${OPT_DIR}/apache-ant-1.9.13
if [[ ! -d ${ANT_DIR} ]]; then
wget -O- https://archive.apache.org/dist/ant/binaries/apache-ant-1.9.13-bin.tar.gz | tar zxvf -
fi

BZIP_DIR=${OPT_DIR}/bzip2-1.0.6
if [[ ! -d ${BZIP_DIR} ]]; then
    wget -O- https://www.sourceware.org/pub/bzip2/bzip2-1.0.6.tar.gz | tar zxvf -
    pushd ${BZIP_DIR}
    make install PREFIX=${BZIP_DIR}_install CFLAGS=" -fPIC"
    popd
fi

# Download ART
ART_DIR=${OPT_DIR}/ART
if [[ ! -d ${ART_DIR} ]]; then
    mkdir -p ${ART_DIR}
    pushd ${ART_DIR}
    wget --no-check-certificate -O- http://www.niehs.nih.gov/research/resources/assets/docs/artbinvanillaicecream031114linux64tgz.tgz | tar xfz -
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
