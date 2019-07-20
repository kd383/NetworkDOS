#!/usr/bin/env python
"""
This is a test function for eig_rand1 module

11/29/2016 - Created
11/29/2016 - Change path import

"""

import sys
from os.path import dirname,realpath
sys.path.insert(0,dirname(realpath(__file__))[:-11])
from pyDOS import *
import scipy.io as sio

def test_plot_cheb(A,n):
	L = matrix_normalize(A)
	c = moments_cheb_dos(L,n,N=50)[0]
	cl = moments_cheb_ldos(L,n,N=50)[0]

	y = plot_cheb((c,))
	yl,idx = plot_cheb_ldos((cl,))
	yh = plot_chebhist((c,))
	yi = plot_chebint((c,))
	yp = plot_chebp((c,))

	idx.reshape([n,-1])
	sio.savemat('test_plt_chb.mat',{'c':c,'cl':cl,'y':y,'yl':yl,'idx':idx,
									'yh':yh,'yi':yi,'yp':yp})


if __name__ == '__main__':
	A = load_graph('erdos02-cc','../data/')
	test_plot_cheb(A,A.shape[0])