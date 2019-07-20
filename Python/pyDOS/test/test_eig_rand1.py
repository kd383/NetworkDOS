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
import numpy.linalg as nl

def test_eig_rand1():
	Acore = nr.randn(10,10);
	Acore = (Acore + Acore.transpose())/2;
	[U,R] = nl.qr(nr.randn(1000,10))
	Afun = lambda x: U.dot(Acore.dot(U.transpose().dot(x)))
	Z = nr.randn(1000,10)
	[V,D] = eig_rand1(Z,Afun)
	sio.savemat('test_eig_rand.mat',{'Acore':Acore,'U':U,'Z':Z,'V':V,'D':D})

if __name__ == '__main__':
	test_eig_rand1()