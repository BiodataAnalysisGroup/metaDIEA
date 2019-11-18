# Dockerized Hisat2-StringTie Pipeline

Execute a containerized version of Hisat2-StringTie Pipeline in an ubuntu-based image. 

## Getting Started

The present Dockerfile creates a 1.08GB image which includes all necessary software and extra packages needed in order to execute the Hisat2-StringTie Pipeline for RNA-seq analysis. 

### Prerequisites

Docker version 19.03 or newer.

### Installing

The bold parts of each command should be changed to reflect the user's setup

Download the Dockerfile in someplace visible from your path, preferably inside an empty folder.

Build the image from scratch.


>$ docker build -t hisat2_stringtie_v.1 **.**


Once the installation has finished you can start a container based on this image.


>$ docker run --name=hisat2Container -it hisat2_stringtie_v.1


This will launch the container in interactive mode. You can make sure everything went well by checking all necessary tools are installed. 


> #/ hisat2

> #/ samtools

> #/ stringtie

> #/ R

> \> library("DESeq2")


## Executing the Pipeline

In order for the Pipeline to be executed correctly, your files should have the following structure.

![Structure](https://user-images.githubusercontent.com/56021536/68419821-a6ebf500-01a3-11ea-90b4-3b242ec9be11.png)

Next, run the following Docker command

>$ docker run --name=tophat2Container -v **//c/Users/Konstantinos/Desktop/Thesis/Hisat2Pipeline**:/app hisat2_stringtie_v.1

The -v flag is used to mount the files included in the path specified inside the container's /app folder.

## Versioning

First version of Containerized Hisat2-StringTie Pipeline

## Authors

* **Konstantinos Koukoutegos** 
* **Fotis Psomopoulos** 
* **Alexandros Dimopoulos** 
* **Panagiotis Moulos** 


## License

This project is licensed under the MIT License.



