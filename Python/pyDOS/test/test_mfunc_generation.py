#!/usr/bin/env python
"""
This is a test function for mfunc_generation module

11/28/2016 - Created
11/29/2016 - Change path import

"""

import sys
from os.path import dirname,realpath
sys.path.insert(0,dirname(realpath(__file__))[:-11])
from pyDOS import *
import scipy.io as sio
import numpy.random as nr

def test_mfunc_generation(A, n, x=0):
	if isinstance(x,int):
		x = nr.randn(n,10)
	Afun1 = mfunc_laplacian(A)
	Afun2 = mfunc_modularity(A)
	Afun3 = mfunc_normalize(A)
	Afun4 = mfunc_slaplacian(A)
	c = nr.randn(10,1)
	Afun5 = mfunc_cheb_poly(c,A)
	sio.savemat('test_mfc_gen.mat',{'c':c, 'x':x,'Afun1x':Afun1(x),'Afun2x':Afun2(x),'Afun3x':Afun3(x),
						'Afun4x':Afun4(x),'Afun5x':Afun5(x)})

if __name__ == '__main__':
	A = load_graph('erdos02-cc','../data/')
	test_mfunc_generation(A, A.shape[0])

