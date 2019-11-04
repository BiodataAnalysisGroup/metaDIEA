# Dockerized Tophat2 Pipeline

Execute a containerized version of Tophat2 Pipeline in an ubuntu-based image. 

## Getting Started

The present Dockerfile creates a 350MB image which includes all necessary software and extra packages needed in order to execute the Tophat2 Pipeline for RNA-seq analysis. 

### Prerequisites

Docker version 19.03 or newer.

### Installing

Download the Dockerfile in someplace visible from your path, preferably inside an empty folder.

Build the image from scratch.


>$ docker build -t tophat2_v.1 .


Once the installation has finished you can start a container based on this image.


>$ docker run --name=tophat2Container -it tophat2_v.1


This will launch the container in interactive mode. You can make sure everything went well by checking all necessary tools are installed. 


> #/ bowtie

> #/ tophat2

> #/ cufflinks


## Executing the Pipeline

In order for the Pipeline to be executed correctly, your files should have the following structure.

![structure](https://user-images.githubusercontent.com/56021536/68129931-23bc6c00-ff23-11e9-94c5-27aa150beb84.png)

Next, run the following Docker command

>$ docker run --name=tophat2Container -v **//c/Users/Konstantinos/Desktop/Thesis/Tophat2Pipeline**:/app tophat2_v.1

The -v flag is used to mount the files included in the path specified inside the container's /app folder.

## Versioning

First version of Containerized Tophat2 Pipeline

## Authors

* **Konstantinos Koukoutegos** 
* **Fotis Psomopoulos** 
* **Alexandros Dimopoulos** 
* **Panagiotis Moulos** 


## License

This project is licensed under the MIT License.
