# Meta DE Analysis

Shiny application for machine learning model building and prediction used for Differential Expression analysis

## Getting Started

Meta DE Analysis is an application written in R programming language, built under the Shiny framework. 

### Prerequisites

R 4.0.0 

### Installing

Once you've downloaded the source code and files make sure first to install all necessary packages used by the application.

Inintially, download the BiocManager, a package manager used by Bioconductor.  

```
$ if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

Once the installation has finished you can start installing all other necessary software.

```
$ install.packages(c("shiny", "randomForest", "xgboost", "ROCR", "MLmetrics", "caret", "magrittr", "dplyr", "Matrix", "pracma"))
```

During the installation you might be asked to update some specific packages. Update all.


## Executing the application

In your R Studio IDE open the app.R file.

On the right top corner of the script window press the "Run App" option. Your console will start printing some information about the libraries being using as well as for the data being processed. 

Please be patient. The application window will launch in a moment. 

## Versioning

Meta DE Analysis 1.0.0

## Authors

* **Konstantinos Koukoutegos** 
* **Fotis Psomopoulos** 
* **Alexandros Dimopoulos** 
* **Panagiotis Moulos** 


## License

This project is licensed under the MIT License.
