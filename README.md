Network Density of States (NDOS)
===============

This repository contains the experiments of the [Network Density of States](https://arxiv.org/abs/1905.09758) by Kun Dong, Austin Benson, David Bindel. 
This paper will be appearing at KDD 2019.

The bibliographic information for the paper will be updated once the proceeding 
is published.
```bibtex
@inproceedings{dong2019networkdos,
  title={Network Density of States},
  author={Dong, Kun and Benson, Austin R Bindel, David},
  booktitle={Proceedings of the 25th ACM SIGKDD International Conference on 
  	Knowledge Discovery \& Data Mining},
  year={2019},
  organization={ACM}
}
```

## Contents
1. [Introduction](#introduction)
2. [Setup](#setup)
3. [Usage](#usage)
	1. [A](#a)
	2. [B](#b)
	3. [C](#b)

## Introduction

Spectral analysis connects graph structure to the eigenvalues and eigenvectors 
of associated matrices.  Much of spectral graph theory has focused on results 
involving only a few extreme eigenvalues and their associated eigenvalues; The 
study of graphs through the overall distribution of eigenvalues --- 
_spectral density_ --- is largely limited to simple random graph models; and 
the interior of the spectrum of real-world graphs remains largely unexplored, 
difficult to compute and to interpret.

In our paper, we borrow tools developed in condensed matter physics, and add 
novel adaptations to handle the spectral signatures of common graph motifs. The 
resulting methods are highly efficient, as we are able to compute spectral 
densities for graphs with over a billion edges even on a single compute node.

## Setup

Clone the repository.

For **Matlab**:
  * If you would like to download the data used in our demos, check the data 
  folder where we supply a shell script _fetch.sh_.  
  * run _startup.m_.

For **Python**, 

## Usage

### A

### B

### C

