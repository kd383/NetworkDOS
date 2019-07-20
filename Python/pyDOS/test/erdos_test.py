#!/usr/bin/env python
import sys
from os.path import dirname,realpath
sys.path.insert(0,dirname(realpath(__file__))[:-11])
from pyDOS import *
from pyDOS.test import *
import numpy as np
import numpy.random as nr
import numpy.linalg as nl
import scipy as sp
import scipy.sparse as ss
import scipy.io as sio
import scipy.sparse.linalg as sl
from copy import deepcopy




if __name__ == '__main__':

	A = load_graph('erdos02-cc','../data/')
	n = A.shape[0]

	test_matrix_generation(A)

	test_mfunc_generation(A,n)

	test_eig_rand1()

	test_moment_filter()

	test_rescale_matrix(A,n)

	test_moment_comp(A,n)

	test_plot_cheb(A,n)

