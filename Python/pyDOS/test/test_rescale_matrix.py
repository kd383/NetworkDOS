#!/usr/bin/env python
"""
This is a test function for rescale_matrix module

11/29/2016 - Created
11/29/2016 - Change path import

"""

import sys
from os.path import dirname,realpath
sys.path.insert(0,dirname(realpath(__file__))[:-11])
from pyDOS import *
import scipy.io as sio
import numpy.random as nr

def test_rescale_matrix(A,n,x=0):
	if isinstance(x,int):
		x = nr.randn(n,1)

	L = mfunc_normalize(A)
	[As, ab] = rescale_matrix(A)
	[Ls, cd] = rescale_mfunc(L,n)
	sio.savemat('test_rsc_mat.mat',{'x':x,'As':As,'ab':ab,'Lx':Ls(x),'cd':cd})

if __name__ == '__main__':
	A = load_graph('erdos02-cc','../data/')
	test_rescale_matrix(A,A.shape[0])