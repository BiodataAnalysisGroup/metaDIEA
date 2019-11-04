# Dockerized Kallisto Pipeline

Execute a containerized version of Kallisto Pipeline in an ubuntu-based image. 

## Getting Started

The present Dockerfile creates a 950MB image which includes all necessary software and extra packages needed in order to execute the Kallisto Pipeline for RNA-seq analysis. 

### Prerequisites

Docker version 19.03 or newer.

### Installing

Download the Dockerfile in someplace visible from your path, preferably inside an empty folder.

Build the image from scratch.

```
$ docker build -t image_name .
```

Once the installation has finished you can start a container based on this image.

```
$ docker run --name=container_name -it image_name
```

This will launch the container in interactive mode. You can make sure everything went well by checking all necessary tools are installed. 

```
#/ kallisto
```

```
#/ R
> library("sleuth")
>q()
```

## Executing the Pipeline

In order for the Pipeline to be executed correctly, your files should have the following structure.

<img src="C:\Users\Konstantinos\Desktop\Thesis\TestKallisto\Structure.png" alt="Structure" title="File Structure" width="221" height="152" />

Next, run the following Docker command

```
$ docker run --name=container_name -v /absolute/path/to/folder/resources:/app image_name
```
The -v flag mounts the files in the path specified inside the container's /app folder.

## Versioning

First version of Containerized Kallisto Pipeline

## Authors

* **Konstantinos Koukoutegos** 
* **Fotis Psomopoulos** 
* **Alexandros Dimopoulos** 
* **Panagiotis Moulos** 


## License

This project is licensed under the MIT License.



