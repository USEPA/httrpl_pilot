# HTTr Pipeline Code (httrpl)

## Installation

### System Requirements

This package provides scripts and APIs in both python 3 and R.  For running the level 0 processing steps (handling of fastq files), python 3 is required.  For running DESeq2, R is required.

#### Python 3

All code in this package has been tested on Python 3.6 - compatibility with newer and older versions has not been tested.

Python 3.6 support can be provided on Red Hat 7 systems by installing the following packages:
```
python36
python36-devel
python36-pip
```

Alternatively, an individual user can install Python 3 on Linux or Windows using [Miniconda](https://docs.conda.io/en/latest/miniconda.html). 

Once Python 3 and associated pip tool is installed, the following python libraries should be installed for the individual user:
```bash
pip3 install --user pandas pymongo mongoengine
```

#### R Install

All code in this package has been tested on R 3.6.0, but any version of R 3.x should be sufficient. 
Package versions may have an impact due to changes in functionality or API.  
The following R libraries are also required to use all features of the R code in Rlib/ directory:

+ [mongolite](https://cran.r-project.org/web/packages/mongolite/index.html) - Required for DB access in R, code has been tested with v2.1.0.
+ [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html) - Required for DB access and handling JSON files, code has been tested with v1.6.
+ [reldist](https://cran.r-project.org/web/packages/reldist/index.html) - Required for computing Gini coefficients when applying QC checks on Level 1 count data, code has been tested with v1.6-6.
+ [DESeq2](https://bioconductor.org/packages/3.9/bioc/html/DESeq2.html) - Required for differential expression analysis, code has been tested and all published analysis run with v1.24.0.


### Installing the httrpl_pilot packages

Currently, it is recommended to clone the entire repository to a user or analysis directory by running:
```bash
git clone https://github.com/USEPA/httrpl_pilot.git (/path/to/analysis)
```

To confirm that all requirements have been successfully installed, run the following command:
```bash
python3 (localpath)/httrpl/bin/testInstall.py
```


## Running the Pipeline

Executable scripts are in: `bin/`
Python modules are in: `lib/`
R modules are in: `Rlib/`


## Contributors

+ **[Logan J. Everett](mailto:everett.logan@epa.gov)**
+ **Imran Shah**
+ **Joseph Bundy**
+ **Derik Haggard**
+ **Beena Vallanat**
