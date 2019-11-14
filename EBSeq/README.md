# Dockerized EBSeq Pipeline

Execute a containerized version of EBSeq Pipeline in an ubuntu-based image. 

## Getting Started

The present Dockerfile creates a 957MB image which includes all necessary software and extra packages needed in order to execute the EBSeq Pipeline for RNA-seq analysis. 

### Prerequisites

Docker version 19.03 or newer.

### Installing

Download the Dockerfile in someplace visible from your path, preferably inside an empty folder.

Build the image from scratch.


>$ docker build -t ebseq_v.1 .


Once the installation has finished you can start a container based on this image.


>$ docker run --name=ebseqContainer -it ebseq_v.1


This will launch the container in interactive mode. You can make sure everything went well by checking all necessary tools are installed. 

> #/ R

> \> library("EBSeq")


## Executing the Pipeline

In order for the Pipeline to be executed correctly, your files should have the following structure.

![structure](https://user-images.githubusercontent.com/56021536/68848700-69262980-06d9-11ea-8d75-72a890669523.png)

Next, run the following Docker command

>$ docker run --name=ebseqContainer -v **//c/Users/Konstantinos/Desktop/Thesis/EBSeqPipeline**:/app ebseq_v.1

The -v flag is used to mount the files included in the path specified inside the container's /app folder.

## Versioning

First version of Containerized EBSeq Pipeline

## Authors

* **Konstantinos Koukoutegos** 
* **Fotis Psomopoulos** 
* **Alexandros Dimopoulos** 
* **Panagiotis Moulos** 


## License

This project is licensed under the MIT License.



