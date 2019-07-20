#!/usr/bin/env python
"""
This module is a function to approximate eigenpairs of approximately low-rank 
symmetric matrix.

8/6/2016 - Created
8/16/2016 - Commented
11/28/2016 - Tested
"""

import numpy as np
import numpy.linalg as npla

def eig_rand1(Z,AZ,thresh=1e-8):
	"""
	If A is an (approximately) low-rank symmetric matrix, estimate eigenpairs 
	of A associated from multiplication by a random probe matrix. Uses a one-pass
	approach -- this could be more accurate with repeated multiplication by A.

	Args:
		Z: A random probe matrix
		AZ: A*Z
		Afun: A function used to apply A to a block of vectors
		thresh: Singular value decay threshold used to truncate A*Z 

	Output:
		V: Eigenvectors sorted by scale of corresponding eigenvalues
		D: A vector of eigenvalues sorted bt scale

	Ref:
		Uses algorithm 5.6 from Halko, Martinsson, and Tropp.
		(http://arxiv.org/pdf/0909.4061v2.pdf)
	"""
	
	# Apply the function
	if callable(AZ):
		AZ = AZ(Z)

	# Get Q an approximate basis for the range space
	(U,s,V) = npla.svd(AZ,0)
	m = np.sum(s>thresh*s[0])
	Q = U[:,np.arange(m)]

	# Should have B*(Q'*z) =Q'*Az; choose B to minimize least square diff
	# The stationary equations look like
	# 		symm(B*M*M'-N*M') = 0
	# or
	# 		B*M*M'+M*M'*B = N*M'+M*N'
	# where M = Q'*Z and  N = Q'*AZ. This is a Sylvester equation, which we
	# solve via a Bartels-Stewart approach

	M = np.dot(Q.T,Z)
	N = np.dot(Q.T,AZ)
	C = np.dot(N,M.T)
	C = C+C.T

	(U,s2,V) = npla.svd(M.T,0)
	V = V.T
	s2 = s2**2
	s2.shape = (len(s2),1)
	e = np.ones(s2.shape)
	Bt = np.dot(np.dot(V.T,C),V)/(np.dot(s2,e.T)+np.dot(e,s2.T))
	B = np.dot(np.dot(V,Bt),V.T)

	# Reconstruct the eigenpairs
	D,V = npla.eig(B)
	I = np.argsort(np.absolute(D))
	I = I[::-1]
	D = D[I]
	V = np.dot(Q,V[:,I])

	return V,D

if __name__ == '__main__':
	pass
