# PrimerTree #
PrimerTree: Visually Assessing the Specificity and Informativeness of Primer Pairs

* [Features](#features)
* [Examples](#examples)
* [Installation](#installation)
* [Usage](#usage)

## Features ##
* Automatically query a primer pair and generate a tree of the products
* Can be run and plotted in two commands
* All options of Primer-Blast and clustal omega are supported.

## Installation ##
### R installation ###
```r
install.packages('primerTree')
```
### Clustal omega installation ###
[Clustal Omega](http://www.clustal.org/omega/#Download)
For windows use the precompiled binary

Simple installation
```bash
./configure && make && make install
```

## Usage ##
Simple search for a Mammal 16S primer
```r
library(primerTree)
mammals_16S = search_primer_pair(name='Mammals 16S', 'CGGTTGGGGTGACCTCGGA', 'GCTGTTATCCCTAGGGTAACT')
plot(mammals_16S)
```

Using the parallel features with the multicore package using the doMC backend, with 8 threads.
```r
library(doMC)
registerDoMC(8)
library(primerTree)
mammals_16S = search_primer_pair(name='Mammals 16S', 'CGGTTGGGGTGACCTCGGA', 'GCTGTTATCCCTAGGGTAACT', .parallel=T)
plot(mammals_16S)
```
