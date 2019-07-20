#!/usr/bin/env python
"""
This module is a colletion of functions that rescale the symmetrix matrix H so 
the eigenvalue range between -1 and 1. 

8/3/2016 - Created
8/16/2016 - Commented
11/28/2016 - Tested
"""

import numpy as np
import scipy.sparse as ss
import scipy.sparse.linalg as ssla

def rescale_matrix(H,range=0.01):
	"""
	Rescale the symmetric matrix H so the eigenvalue range maps to between -1 
	and 1. If range is not given, the function uses Lanczos to estimate the 
	extremal eigenvalues, expanding by a relative fudge factor.

	Args:
		H: The original matrix
		range: A two-element vector representing an interval of eigenvalues
		fudge: A scalar to guard against range underestimates (default=0.01)
	
	Output:
		H: The scaled matrix
		ab: Transformation parameters: Hs=(H-b)/a
	"""

	# Get the fudge factor and range
	if not isinstance(range,list) or len(range) == 1:
		fudge = range
		range = []
	else:
		fudge = 0

	n = H.shape[0]

	# Compute range if not given
	if not range:
		range = ssla.eigsh(H,1,which='SA',
				return_eigenvectors=False)
		range = np.append(range,ssla.eigsh(H,1,
				which='LA',return_eigenvectors=False))

	# Form dense or sparse identity, as appropraite
	if ss.issparse(H):
		I = ss.eye(H.shape[0])
	else:
		I = np.eye(H.shape[0])

	# Parameters for linear mapping
	ab = [(range[1]-range[0])/(2-fudge),(range[1]+range[0])/2]
	H = (H-ab[1]*I)/ab[0]
	ab = np.asarray(ab).reshape([2,-1])
	return H,ab

def rescale_mfunc(H,n=0.01,range=0.01):
	"""
	Rescale the symmetric matrix/operator H so the eigenvalue range maps to 
	between -1 and 1. If the range is not given, the function uses Lanczos
	to estimate the extremal eigenvalues, expanding by a relative fudge factor

	Input:
		H: The original operator/matrix
		n: The dimension of the space
		range: A two-element vector representing an interval of eigenvalues
		fudge: A scalar to guard against range underestimates (default=0.01) 

	Output:
		Hfun: The scaled operator
		ab: Transformation parameters: Hs=(H-b)/a
	"""
	
	# Create a function handle if a matrix is given
	if callable(H):
		if not isinstance(n,int):
			raise ValueError('Missing size argument')
		if not isinstance(range,list) or len(range) == 1:
			fudge = range
			range = []
		else:
			fudge = 0

		# Lanczos to estimate range if not given
		if not range:
			range = ssla.eigsh(ssla.LinearOperator((n,n),H,dtype='float64'),1,which='SA',
				return_eigenvectors=False)
			range = np.append(range,ssla.eigsh(ssla.LinearOperator((n,n),H,dtype='float64'),1,
				which='LA',return_eigenvectors=False))

		# Parameters for lienar mapping
		ab = [(range[1]-range[0])/(2-fudge),(range[1]+range[0])/2]
		Hfun = lambda x: (H(x)-ab[1]*x)/ab[0]
	else:
		if isinstance(H,np.ndarray):
			H = ss.csr_matrix(H)
		if not isinstance(n,int):
			range = n
			n = H.shape[0]
		H,ab = rescale_matrix(H,range)
		Hfun = lambda x: H*x

	ab = np.asarray(ab).reshape([2,-1])
	return Hfun,ab

if __name__ == '__main__':
	pass
