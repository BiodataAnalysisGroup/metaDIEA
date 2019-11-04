FROM ubuntu:latest

MAINTAINER konstantinos.koukoutegos@gmail.com

#Update, download Kallisto binary and copy to /bin/
RUN	apt-get update \
	&& apt install wget -y \
	&& wget https://github.com/pachterlab/kallisto/releases/download/v0.46.0/kallisto_linux-v0.46.0.tar.gz \ 
	&& tar xvzf kallisto_linux-v0.46.0.tar.gz \ 
	&& rm kallisto_linux-v0.46.0.tar.gz \
	&& cp -r /kallisto/kallisto /usr/local/bin/ \ 
	&& rm -r /kallisto
	
# Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

ENV DEBIAN_FRONTEND noninteractive

#Update and install r-base
RUN apt-get update \ 
	&& apt-get install -y --no-install-recommends \
	    apt-utils \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		apt-transport-https \
		gsfonts \
		gnupg2 \
		r-base \ 
	&& rm -rf /var/lib/apt/lists/*

#install devtools and sleuth
RUN apt-get update && \
    apt-get install -y libcurl4-openssl-dev libssl-dev libssh2-1-dev libxml2-dev && \
    R -e "install.packages(c('devtools', 'testthat', 'roxygen2'))" && \
	  R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')" && \
	  R -e "BiocManager::install('rhdf5')" && \
	  R -e "devtools::install_github('pachterlab/sleuth')" 
	
#run kallisto.log.txt
CMD ["/bin/bash", "app/readme.kallisto.log.txt"]	

