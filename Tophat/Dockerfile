FROM ubuntu:latest

MAINTAINER konstantinos.koukoutegos@gmail.com

#install bowtie
RUN apt-get update \
    && apt-get install bowtie -y \
	  && apt install wget -y 
	
#install tophat-2
RUN wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz \
	&& tar xvzf tophat-2.1.1.Linux_x86_64.tar.gz \
	&& rm tophat-2.1.1.Linux_x86_64.tar.gz \
	&& cp -r /tophat-2.1.1.Linux_x86_64/. usr/local/bin \
	&& rm -r /tophat-2.1.1.Linux_x86_64
	
#install cufflinks
RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz \
	&& tar xvzf cufflinks-2.2.1.Linux_x86_64.tar.gz \
	&& rm cufflinks-2.2.1.Linux_x86_64.tar.gz \
	&& cp -r /cufflinks-2.2.1.Linux_x86_64/. usr/local/bin \ 
	&& rm -r cufflinks-2.2.1.Linux_x86_64 


#run tophat2_cufflinks.log.txt
CMD ["/bin/bash", "app/readme.tophat2_cufflinks.log.txt"]
