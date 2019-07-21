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
	1. [Compute DOS/PDOS](#compute-dospdos)
	2. [Motif Filtering](#motif-filtering)
	3. [Model Comparison: BTER](#model-comparison-bter)
	4. [Other Demos](#other-demos)

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

* `pyDOS` is implemented in Python2. Demos are coming soon.

## Usage

### Compute DOS/PDOS

Our first demo is a simple computation of density of states (DOS) and pointwise
density of states (PDOS). To run this experiment and produce the figure below, 
you can use the `demo_dos` and `demo_ldos` commands. By default, they use 
**Kernel Polynomial Method (KPM)** to compute 1000 Chebyshev moments with 20 
probe vectors for the Erd&#337;s Collaboration Network.

* method: 'cheb'(default), 'lan', 'nd', 'exact'
* dname: Use one from RODGER dataset, or supply an adjacency matrix.

Note that 'exact' method uses full matrix-matrix multiplication, so it should 
be avoided on large networks. When 'lan' is used for PDOS, we only compute a 
subsample of 100 nodes.

<p align="center">
    <img src="/pics/erdos_dos.png" width="220">
    <img src="/pics/erdos_dos_zoom.png" width="220">
    <img src="/pics/erdos_ldos.png" width="220">
</p>

For DOS, the blue bars are the exact count of eigenvalues in each bin, and the
red dots are our approximation. The middle figure zooms in near the bottom. For 
PDOS, the y-axis represents the node index, the x-axis represents eigenvalues, 
and the colors indicate the heights of the spectral histogram.

### Motif Filtering
<p align="center">
    <img src="/pics/hepth_0filt.png" width="300">
    <img src="/pics/hepth_1filt.png" width="300">
    <br>
    (a) No Filter &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; (b) Filter at &#955;=0
    <br>
    <img src="/pics/hepth_2filt.png" width="300">
    <img src="/pics/hepth_3filt.png" width="300">
    <br>
    (c) Filter at &#955;=-1/3 &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; (d) Filter at &#955;=-1/2
    <br>
    <img src="/pics/hepth_4filt.png" width="300">
    <img src="/pics/filter_error.png" width="300">
    <br>
    (c) Filter at &#955;=-1/3 &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; (d) Relative Error
</p>

### Model Comparison: BTER

[BTER model](https://arxiv.org/pdf/1112.3644.pdf) closely captures many properties of a graph, such as degree distribution, clustering coefficient, and distribution of eigenvalues of the adjacency matrix. We followed this [guide](https://www.sandia.gov/~tgkolda/feastpack/doc_bter_match.html) to create a model for the Erd&#337;s collaboration network and compare its DOS/PDOS against the original. We find the construction process of BTER, particular its treatment of degree-one nodes, leads to an abundance of motifs absent in the original graph. Use the command `demo_bter_comparison` to run this demo.

<p align="center">
    <img src="/pics/erdos_dos_zoom.png" width="300">
    <img src="/pics/bter_dos.png" width="300">
    <br>
    (a) Erd&#337;s DOS &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; (b) BTER DOS
    <br>
    <img src="/pics/erdos_ldos.png" width="300">
    <img src="/pics/bter_ldos.png" width="300">
    <br>
    (c) Erd&#337;s PDOS &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; (d) BTER PDOS
</p>

### Other demos

A few demos are coming in the near future. Thank you for your patience.