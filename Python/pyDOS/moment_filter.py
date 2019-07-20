#!/usr/bin/env python
"""
This module is a colletion of functions that apply filters on Chebyshev moments

8/8/2016 - Created
8/16/2016 - Commented
11/28/2016 - Tested
"""

import numpy as np
import scipy.sparse as ss

def filter_jackson(c):
	"""
	Apply the Jackson filter to a sequence of Chebyshev	moments. The moments 
	should be arranged column by column.

	Args:
		c: Unfiltered Chebyshev moments

	Output:
		cf: Jackson filtered Chebyshev moments
	"""

	N = c.shape[0]
	n = np.arange(N)
	tau = np.pi/(N+1)
	g = ((N-n+1)*np.cos(tau*n)+np.sin(tau*n)/np.tan(tau))/(N+1)
	g.shape = (N,1)
	c = g*c
	return c

def filter_lorentz(c,beta=4.0):
	"""
	Apply the Lorentz filter to a sequence of Chebyshev	moments. The moments 
	should be arranged column by column.

	Args:
		c: Unfiltered Chebyshev moments

	Output:
		cf: Lorentz filtered Chebyshev moments
	"""

	N = c.shape[0]
	n = np.arange(N,dtype='float')
	g = np.sinh(beta*(1-n/N))/np.sinh(beta)
	g.shape = (N,1)
	c = g*c
	return c

if __name__ == '__main__':
	pass