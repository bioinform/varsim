FROM centos:7.2.1511
		
ARG branch_of_interest		
ENV VARSIM_VERSION $branch_of_interest

RUN yum update -y
RUN yum install -y centos-release-scl
RUN yum install -y devtoolset-4
RUN yum install -y build-essential zlib1g-dev curl less vim bzip2
RUN yum install -y git wget

RUN cd /opt && git clone https://github.com/bioinform/varsim.git
RUN cd /opt/varsim && git checkout ${VARSIM_VERSION} && ./build.sh
ENV PATH=/opt/varsim/:${PATH}
