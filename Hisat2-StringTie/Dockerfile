FROM ubuntu:latest

MAINTAINER konstantinos.koukoutegos@gmail.com

#install hisat2, samtools and wget
RUN apt-get update \
	&& apt-get install hisat2 -y \
	&& apt install samtools -y \
	&& apt install wget -y 
	
#install stringtie-2	
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.0.4.Linux_x86_64.tar.gz \
	&& tar xvzf stringtie-2.0.4.Linux_x86_64.tar.gz \
	&& rm stringtie-2.0.4.Linux_x86_64.tar.gz \
	&& cp -r /stringtie-2.0.4.Linux_x86_64/. usr/local/bin \
	&& rm -r /stringtie-2.0.4.Linux_x86_64
	
# Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

ENV DEBIAN_FRONTEND noninteractive

#install r 3.6.1
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" >> etc/apt/sources.list \
	#&& apt-get update \
	&& apt-get install gnupg2 -y \
	&& apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
	&& apt-get update \
	&& apt-get upgrade -y \
	&& apt-get install r-base-dev -y \
	&& apt-get install libxml2-dev -y \
	&& apt-get install libcurl4-gnutls-dev -y 

#install DESeq2
RUN	R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')" && \
	R -e "BiocManager::install('DESeq2')" 

#run scripts
CMD ["/bin/bash", "app/genome/make_mm10.sh"]
