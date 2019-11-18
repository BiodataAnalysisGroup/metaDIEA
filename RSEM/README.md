# Dockerized RSEM Pipeline

Execute a containerized version of RSEM Pipeline in an ubuntu-based image. 

## Getting Started

The present Dockerfile creates a 808MB image which includes all necessary software and extra packages needed in order to execute the RSEM Pipeline for RNA-seq analysis. 

### Prerequisites

Docker version 19.03 or newer.

### Installing

The bold parts of each command should be changed to reflect the user's setup.

Download the Dockerfile in someplace visible from your path, preferably inside an empty folder.

Build the image from scratch.


>$ docker build -t rsem_v.1 **.**


Once the installation has finished you can start a container based on this image.


>$ docker run --name=rsemContainer -it rsem_v.1


This will launch the container in interactive mode. You can make sure everything went well by checking all necessary tools are installed. 

> #/ rsem-prepare-reference

> #/ rsem-run-ebseq --ngvector

> #/ bowtie

> #/ samtools

> #/ R

> \> library("EBSeq")


## Executing the Pipeline

In order for the Pipeline to be executed correctly, your files should have the following structure.

![Structure](https://user-images.githubusercontent.com/56021536/68848480-0d5ba080-06d9-11ea-9a46-e6609ef9dfc1.png)

Next, run the following Docker command

>$ docker run --name=rsemContainer -v **//c/Users/Konstantinos/Desktop/Thesis/RSEMPipeline**:/app rsem_v.1

The -v flag is used to mount the files included in the path specified inside the container's /app folder.

## Versioning

First version of Containerized RSEM Pipeline

## Authors

* **Konstantinos Koukoutegos** 
* **Fotis Psomopoulos** 
* **Alexandros Dimopoulos** 
* **Panagiotis Moulos** 


## License

This project is licensed under the MIT License.



