#!/usr/bin/env python
"""
This module is a colletion of functions that compute the Chebyshev moments 
(density of states) of various operators. 

8/6/2016 - Created
8/16/2016 - Commented
"""

import numpy as np
import scipy.sparse as ss
import scipy.sparse.linalg as ssla
import numpy.random as nr

def moments_cheb(A,V,N=10,kind=1):
	"""
	Compute a column vector of Chebyshev moments of the form c(k) = v'*T_k(A)*v 
	for k = 0 to N-1. This routine does no scaling; the spectrum of A should 
	already lie in [-1,1]

	Args:
		A: Matrix or function apply matrix (to multiple RHS)
		V: Starting vectors
		N: Number of moments to compute
		kind: 1 or 2 for first or second kind Chebyshev functions
			(default = 1)

	Output:
		c: a length N vector of moments
	"""

	if N<2:
		N = 2

	if not isinstance(V,np.ndarray):
		V = V.toarray()

	# Create a function handle if given a matrix
	if  callable(A):
		Afun = A
	else:
		if isinstance(A,np.ndarray):
			A = ss.csr_matrix(A)
		Afun = lambda x: A*x

	n,p = V.shape
	c = np.zeros((N,p))

	# Run three-term recurrence to compute moments
	TVp = V
	TVk = kind*Afun(V)
	c[0] = np.sum(V*TVp,0)
	c[1] = np.sum(V*TVk,0)
	for i in range(2,N):
		TV = 2*Afun(TVk) - TVp
		TVp = TVk
		TVk = TV
		c[i] = sum(V*TVk,0)

	return c

def moments_cheb_dos(A,n,nZ=100,N=10,kind=1):
	"""
	Compute a column vector of Chebyshev moments of the form c(k) = tr(T_k(A)) 
	for k = 0 to N-1. This routine does no scaling; the spectrum of A should 
	already lie in [-1,1]. The traces are computed via a stochastic estimator 
	with nZ probe

	Args:
		A: Matrix or function apply matrix (to multiple RHS)
		n: Dimension of the space
		nZ: Number of probe vectors with which we compute moments
		N: Number of moments to compute
		kind: 1 or 2 for first or second kind Chebyshev functions
		 	(default = 1)

	Output:
		c: a column vector of N moment estimates
		cs: standard deviation of the moment estimator 
			(std/sqrt(nZ))
	"""

	# Create a function handle if given a matrix 
	if callable(A):
		Afun = A
	else:
		if isinstance(A,np.ndarray):
			A = ss.csr_matrix(A)
		Afun = lambda x: A*x

	if N < 2:
		N = 2

	# Set up random probe vectors (allowed to be passed in)
	if not isinstance(nZ,int):
		Z = nZ
		nZ = Z.shape[1]
	else:
		Z = np.sign(nr.randn(n,nZ))

	# Estimate moments for each probe vector
	cZ = moments_cheb(Afun,Z,N,kind)
	c = np.mean(cZ,1)
	cs = np.std(cZ,1,ddof=1)/np.sqrt(nZ)

	c = c.reshape([N,-1])
	cs = cs.reshape([N,-1])
	return c,cs

def moments_cheb_ldos(A,n,nZ=100,N=10,kind=1):
	"""
	Compute a column vector of Chebyshev moments of the form 
	c(k,j) = [T_k(A)]_jj for k = 0 to N-1. This routine does no scaling; the 
	spectrum of A should already lie in [-1,1]. The diagonal entries are 
	computed by a stochastic estimator

	Args:
		A: Matrix or function apply matrix (to multiple RHS)
		n: Dimension of the space
		nZ: Number of probe vectors with which we compute moments
		N: Number of moments to compute
		kind: 1 or 2 for first or second kind Chebyshev functions
		 	(default = 1)

	Output:
		c: a (N,n) matrix of moments
		cs: standard deviation of the moment estimator 
			(std/sqrt(nZ))
	"""
	
	# Create a function handle if given a matrix
	if callable(A):
		Afun = A
	else:
		if isinstance(A,np.ndarray):
			A = ss.csr_matrix(A)
		Afun = lambda x: A*x

	if N < 2:
		N = 2

	# Set up random probe vectors (allowed to be passed in)
	if not isinstance(nZ,int):
		Z = nZ
		nZ = Z.shape[1]
	else:
		Z = np.sign(nr.randn(n,nZ))

	# Run three-term recurrence to estimate moments.
	# Use the stochastic diagonal estimator of Bekas and Saad
	# http://www-users.cs.umn.edu/~saad/PDF/usmi-2005-082.pdf

	c = np.zeros((N,n))
	cs = np.zeros((N,n))

	TZp = Z
	X = Z*TZp
	c[0] = np.mean(X,1).T
	cs[0] = np.std(X,1,ddof=1).T

	TZk = kind*Afun(Z)
	X = Z*TZk
	c[1] = np.mean(X,1).T
	cs[1] = np.std(X,1,ddof=1).T

	for i in range(2,N):
		TZ = 2*Afun(TZk) - TZp
		TZp = TZk
		TZk = TZ
		X = Z*TZk
		c[i] = np.mean(X,1).T
		cs[i] = np.std(X,1,ddof=1).T

	cs = cs/np.sqrt(nZ)

	c = c.reshape([N,-1])
	cs = cs.reshape([N,-1])
	return c,cs

def moments_delta(x,N):
	"""
	Compute Chebyshev moment 0 through N-1 of a delta distribution centered 
	at x.

	Args:
		x: Center of delta distribution
		N: Number of moments to compute

	Output:
		c: a column vector of N Chebyshev moments
	"""

	c = np.sum(np.real(np.cos(np.linspace(0,N-1,N)*np.arccos(np.complex_(x)))),0)
	c = c.reshape([N,-1])
	return c

if __name__ == '__main__':
	pass