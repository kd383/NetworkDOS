#!/usr/bin/env python
"""
This module is a colletion of functions that generates different matrices for graph

8/1/2016 - Created
8/16/2016 - Commented
11/25/2016 - Tested
11/29/2016 - Fix bug: deepcopy
"""

import numpy as np
import copy as cp
import scipy.sparse as ss

def matrix_adjacency(W,keep_diag=1):
	"""
	Convert a weighted adjacency matrix into an ordinary adjacency.

	Args:
		W: weighted adjacency matrix
		keep_diag: flag if we want to keep the diagonal on output (default: 1)

	Output: 
		A: a 0-1 adjacency matrix (in sparse form)
	"""

	if isinstance(W,np.ndarray):
		W = ss.csr_matrix(W)

	# Zero out the diagonal if it is not wanted
	if not keep_diag:
		z = np.zeros(W.shape[0])
		W.setdiag(z,0)

	# Reconstruct with only off-diagonals
	[i,j,wij] = ss.find(W)
	A = ss.csr_matrix((np.ones(wij.size),(i,j)),shape=W.shape)
	return A


def matrix_laplacian(W,mode='r'):
	"""
	Convert a weighted adjacency matrix into a Laplacian

	Args:
		W: weighted adjacency matrix
		mode: 'r'ow or 'c'ol sum (default is rows)

	Output:
		L: a graph Laplacian
	"""

	if isinstance(W,np.ndarray):
		W = ss.csr_matrix(W)

	# Strip diagonal
	z = np.zeros(W.shape[0])
	L = cp.deepcopy(W)
	L.setdiag(z,0)

	# Compute row or column sums
	if mode == 'r':
		d = np.asarray(L.sum(1))
	else:
		d = np.asarray(L.sum(0))

	# Replace the main diagonal, flip off-diagonal signs
	L = -L
	L.setdiag(d.squeeze(),0)
	return L

def matrix_normalize(W,mode='s'):
	"""
	Normalize a weighted adjacency matrix.

	Args:
		W: weighted adjacency matrix
		mode: string indicating the style of normalization;
			's': Symmetric scaling by the degree (default)
			'r': Normalize to row-stochastic
			'c': Normalize to col-stochastic

	Output:
		N: a normalized adjacency matrix or stochastic matrix (in sparse form)
	"""

	dc = np.asarray(W.sum(0)).squeeze()
	dr = np.asarray(W.sum(1)).squeeze()
	[i,j,wij] = ss.find(W)

	# Normalize in desired style
	if mode in 'sl':
		wij = wij/np.sqrt(dr[i]*dc[j])
	elif mode == 'r':
		wij = wij/dr[i]
	elif mode == 'c':
		wij = wij/dc[j]
	else:
		raise ValueError('Unknown mode!')

	N = ss.csr_matrix((wij,(i,j)),shape=W.shape)
	return N

def matrix_slaplacian(W,mode='r'):
	"""
	Convert a weighted adjacency matrix into a signless Laplacian.

	Args:
		W: weighted adjacency matrix
		mode: 'r'ow or 'c'ol sum to zero (default is rows)

	Output:
		L: a signless graph Laplacian
	"""

	if isinstance(W,np.ndarray):
		W = ss.csr_matrix(W)

	# Strip diagonal
	z = np.zeros(W.shape[0])
	L = cp.deepcopy(W)
	L.setdiag(z,0)

	# Compute row or column sums
	if mode == 'r':
		d = np.asarray(L.sum(1))
	else:
		d = np.asarray(L.sum(0))

	# Replace main diagonal
	
	L.setdiag(d.squeeze(),0)
	return L	


if __name__ == '__main__':
	pass


