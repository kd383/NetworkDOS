#!/usr/bin/env python
"""
This is a test function for moment_filter module

11/29/2016 - Created
11/29/2016 -Change path import

"""

import sys
from os.path import dirname,realpath
sys.path.insert(0,dirname(realpath(__file__))[:-11])
from pyDOS import *
import scipy.io as sio
import numpy.random as nr

def test_moment_filter(c=0):
	if isinstance(c,int):
		c = nr.randn(1000,10)

	cj = filter_jackson(c)
	cl = filter_lorentz(c)
	sio.savemat('test_mmt_fil.mat',{'c':c,'cj':cj,'cl':cl})

if __name__ == '__main__':
	test_moment_filter()