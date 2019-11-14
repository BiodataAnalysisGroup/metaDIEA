FROM ubuntu:latest

MAINTAINER konstantinos.konstantinos@gmail.com

#install bowtie2 and wget
RUN apt-get update \
    && apt-get install bowtie2 -y \
	&& apt install samtools -y \
	&& apt install wget -y 
	
# Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

ENV DEBIAN_FRONTEND noninteractive

#install r 3.6.1
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" >> etc/apt/sources.list \
	&& apt-get install gnupg2 -y \
	&& apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
	&& apt-get update \
	&& apt-get upgrade -y \
	&& apt-get install r-base-dev -y \
	&& apt-get install libxml2-dev -y \
	&& apt-get install libcurl4-gnutls-dev -y
	
#install EBSeq 
RUN	R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')" && \
	R -e "BiocManager::install('EBSeq')"

##download and install RSEM source code	
RUN wget https://github.com/deweylab/RSEM/archive/v1.3.1.tar.gz \
	&& tar xvfz v1.3.1.tar.gz \
	&& rm v1.3.1.tar.gz \
	&& cd RSEM-1.3.1 \
	&& make \
	&& cp rsem-* /usr/local/bin/ \
	&& cp WHAT_IS_NEW /usr/local/bin/ \
	&& cp rsem_perl_utils.pm /usr/local/bin/ \
	&& cd EBSeq \
	&& make \
	&& cp rsem-* /usr/local/bin/ \
	&& cd .. \
	&& cd .. \
	&& rm -r RSEM-1.3.1
	
#run readme.RSEM.log.txt
CMD ["/bin/bash", "app/readme.RSEM.log.txt"]