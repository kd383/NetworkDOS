#!/usr/bin/env python
"""
This is a test function for matrix_generation module

11/28/2016 - Created
11/29/2016 - Change path import

"""

import sys
from os.path import dirname,realpath
sys.path.insert(0,dirname(realpath(__file__))[:-11])
from pyDOS import *
import scipy.io as sio



def test_matrix_generation(A):
	A1 = matrix_adjacency(A)
	A2 = matrix_laplacian(A)
	A3 = matrix_normalize(A)
	A4 = matrix_slaplacian(A)
	sio.savemat('test_mat_gen.mat',{'A1':A1,'A2':A2,'A3':A3,'A4':A4})

if __name__ == '__main__':
	A = load_graph('erdos02-cc','../data/')
	test_matrix_generation(A)

