# PrimerTree #
[![Travis-CI Build Status](https://travis-ci.org/jimhester/primerTree.svg?branch=master)](https://travis-ci.org/jimhester/primerTree)
PrimerTree: Visually Assessing the Specificity and Informativeness of Primer Pairs

* [Features](#features)
* [Examples](#examples)
* [Installation](#installation)
* [Usage](#usage)

## Features ##
* Automatically query a primer pair and generate a tree of the products.
* Analysis can be run and then plotted in two commands.
* All options of Primer-Blast and clustal omega are supported.

## Installation ##
### R installation ###
#### CRAN ####
```s
install.packages('primerTree')
```
#### Github ####
```s
# install.packages("devtools")
devtools::install_github("jimhester/primerTree")
```
### [Clustal Omega](http://www.clustal.org/omega/#Download) Installation ###
#### Windows ####
Use the pre-compiled windows binary.  Either put the installed clustalo.exe in your path, or pass the path to the executable in the clustal_options option
```s
library(primerTree)
mammals_16S = search_primer_pair(name='Mammals 16S',
  'CGGTTGGGGTGACCTCGGA', 'GCTGTTATCCCTAGGGTAACT', clustal_options=c(exec='C:\Program Files\Clustal Omega\clustalo.exe'))
```
#### Linux ####
Simple installation from source
```bash
./configure && make && make install
```
If the resulting clustalo program is in your path it should be automatically found,
otherwise see the windows instructions on how to specify the path to the
executable.

#### Mac OSX ####
Clustal Omega can be installed from source but requires Xcode CLT and argtable, which are not pre-installed on OSX.

First install Xcode CLT.
```bash
xcode-select --install
```
then install argtable
```bash
brew install argtable
```
Now you can install Clustal Omega from source.
From the clustalo directory:

```bash
./configure && sudo make install
```

## Usage ##
Simple search for a Mammal 16S primer
```s
library(primerTree)
mammals_16S = search_primer_pair(name='Mammals 16S', 'CGGTTGGGGTGACCTCGGA', 'GCTGTTATCCCTAGGGTAACT')
plot(mammals_16S)
```

Using the parallel features with the multicore package using the doMC backend, with 8 threads.
```s
library(doMC)
registerDoMC(8)
library(primerTree)
mammals_16S = search_primer_pair(name='Mammals 16S',
  'CGGTTGGGGTGACCTCGGA', 'GCTGTTATCCCTAGGGTAACT', .parallel=T)
plot(mammals_16S)
```
