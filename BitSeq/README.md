# Dockerized BitSeq Pipeline

Execute a containerized version of BitSeq Pipeline in an ubuntu-based image. 

## Getting Started

The present Dockerfile creates a 981MB image which includes all necessary software and extra packages needed in order to execute the BitSeq Pipeline for RNA-seq analysis. 

### Prerequisites

Docker version 19.03 or newer.

### Installing

The bold parts of each command should be changed to reflect the user's setup.

Download the Dockerfile in someplace visible from your path, preferably inside an empty folder.

Build the image from scratch.


>$ docker build -t bitseq_v.1 **.**


Once the installation has finished you can start a container based on this image.


>$ docker run --name=bitseqContainer -it bitseq_v.1


This will launch the container in interactive mode. You can make sure everything went well by checking all necessary tools are installed. 


> #/ bowtie

> #/ samtools

> #/ R

> \> library("BitSeq")


## Executing the Pipeline

In order for the Pipeline to be executed correctly, your files should have the following structure.

![Structure](https://user-images.githubusercontent.com/56021536/68848149-7abb0180-06d8-11ea-81e0-7567286a8239.png)

Next, run the following Docker command

>$ docker run --name=bitseqContainer -v **//c/Users/Konstantinos/Desktop/Thesis/BitSeqPipeline**:/app bitseq_v.1

The -v flag is used to mount the files included in the path specified inside the container's /app folder.

## Versioning

First version of Containerized BitSeq Pipeline

## Authors

* **Konstantinos Koukoutegos** 
* **Fotis Psomopoulos** 
* **Alexandros Dimopoulos** 
* **Panagiotis Moulos** 


## License

This project is licensed under the MIT License.



