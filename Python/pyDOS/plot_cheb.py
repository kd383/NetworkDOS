#!/usr/bin/env python
"""
This module is a collection of functions used to make plots associated with 
Chebyshev moments

8/9/2016 - Created
8/16/2016 - Commented
"""

import numpy as np
import scipy.sparse.linalg as ssla
import matplotlib.pyplot as plt


def plot_cheb_argparse(npts,c,xx0=-1,ab=np.array([1,0])):
	"""
	Handle argument parsing for plotting routines. Should not be called directly
	by users.

	Args:
		npts: Number of points in a default mesh
		c: Vector of moments
		xx0: Input sampling mesh (original coordinates)
		ab: Scaling map parameters

	Output:
		c: Vector of moments
		xx: Input sampling mesh ([-1,1] coordinates)
		xx0: Input sampling mesh (original coordinates)
		ab: Scaling map parameters
	"""

	if isinstance(xx0,int):
		# only c is given
		xx0 = np.linspace(-1+1e-8,1-1e-8,npts)
		xx = xx0
	else:
		if len(xx0)==2:
			# parameters are c, ab
			ab = xx0
			xx = np.linspace(-1+1e-8,1-1e-8,npts)
			xx0 = ab[0]*xx+ab[1]
		else:
			# parameteres are c, xx0
			xx=xx0

	# All parameters specified
	if not (ab==[1,0]).all():
		xx = (xx0-ab[1])/ab[0]

	return c,xx,xx0,ab

def plot_cheb(varargin,pflag=True):
	"""
	Given a set of first-kind Chebyshev moments, compute the associated density.
	Output a plot of the density function by default

	Args:
		c: Vector of Chebyshev moments (on [-1,1])
		xx: Evaluation points (defaults to mesh of 1001 pts)
		ab: Mapping parameters (default to identity)
		pflag: Option to output the plot

	Output:
		yy: Density evaluated at xx mesh
	"""

	# Parse arguments
	c,xx,xx0,ab = plot_cheb_argparse(1001,*varargin)

	# Run the recurrence
	kind = 1
	N = len(c)
	P0 = xx*0+1
	P1 = kind*xx
	yy =c[0]/(3-kind)+c[1]*xx

	for idx in np.arange(2,N):
		Pn = 2*(xx*P1)-P0
		yy += c[idx]*Pn
		P0 = P1
		P1 = Pn

	# Normalization
	if kind == 1:
		yy = (2/np.pi/ab[0])*(yy/(1e-12+np.sqrt(1-xx**2)))
	else:
		yy = (2/np.pi/ab[0])*(yy*np.sqrt(1-xx**2))

	# Plot by default
	if pflag:
		plt.plot(xx0,yy)
		plt.ion()
		plt.show()
		plt.pause(1)
		plt.clf()
		# input('Press [enter] to continue.')

	yy.reshape([1,-1])
	return yy

def plot_cheb_ldos(varargin,pflag=True):
	"""
	Given a set of first-kind Chebyshev moments, compute the associated local 
	density. Output a plot of the local density functions by default.

	Args:
		c: Vector of Chebyshev moments (on [-1,1])
		xx: Evaluation points (defaults to mesh of 1001 pts)
		ab: Mapping parameters (default to identity)
		pflag: Option to output the plot

	Output:
		yy: Density evaluated at xx mesh (size nnodes-by-nmesh)
		index: Index for spectral re-ordering
	"""

	# Parse arguments
	c,xx,xx0,ab = plot_cheb_argparse(51,*varargin)

	# Run the recurrence to compute CDF
	nmoment,nnodes = c.shape
	txx = np.arccos(xx)
	yy = c[0].reshape((nnodes,1))*(txx-np.pi)/2
	for idx in np.arange(1,nmoment):
		yy = yy +c[idx].reshape((nnodes,1))*np.sin(idx*txx)/idx

	# Difference the CDF to compute histogram
	yy *= -2/np.pi
	yy = yy[:,1:]-yy[:,0:-1]

	# Compute the sorted histogram
	(u,s,v) = ssla.svds(yy,1)
	index = np.argsort(np.squeeze(u))

	# Plot by default
	if pflag:
		fig,ax = plt.subplots()
		yr = np.array([1,nnodes])
		xr = np.array([xx0[0]+xx0[1],xx0[-1]+xx0[-2]],dtype='float')/2
		im = ax.imshow(yy[index,:],extent=np.append(xr,yr),aspect=1.5/nnodes)
		fig.colorbar(im)
		plt.ion()
		plt.show()
		plt.pause(1)
		plt.clf()

	index = np.float64(index.reshape([-1,1]))
	return yy,index

def plot_chebhist(varargin,pflag=True):
	"""
	Given a (filtered) set of first-kind Chebyshev moments, compute the integral
	of the density:
		int_0^s (2/pi)*sqrt(1-x^2)*( c(0)/2+sum_{n=1}^{N-1}c_nT_n(x) )
	Output a histogram of cumulative density function by default.

	Args:
		c: Vector of Chebyshev moments (on [-1,1])
		xx: Evaluation points (defaults to mesh of 21 pts)
		ab: Mapping parameters (default to identity)
		pflag: Option to output the plot

	Output:
		yy: Estimated counts on buckets between xx points
	"""

	# Parse arguments
	c,xx,xx0,ab = plot_cheb_argparse(21,*varargin)

	# Compute CDF and bin the difference
	yy = plot_chebint((c,xx0,ab),pflag=False)
	yy = yy[1:]-yy[:-1]
	xm = (xx0[1:]+xx0[:-1])/2
	
	# Plot by default
	if pflag:
		plt.bar(xm,yy,align='center',width=0.1)
		plt.ion()
		plt.show()
		plt.pause(1)
		plt.clf()

	return yy

def plot_chebint(varargin,pflag=True):
	"""
	Given a (filtered) set of first-kind Chebyshev moments, compute the integral
	of the density:
		int_0^s (2/pi)*sqrt(1-x^2)*( c(0)/2+sum_{n=1}^{N-1}c_nT_n(x) )
	Output a plot of cumulative density function by default.

	Args:
		c: Array of Chebyshev moments (on [-1,1])
		xx: Evaluation points (defaults to mesh of 1001 pts)
		ab: Mapping parameters (default to identity)
		pflag: Option to output the plot

	Output:
		yy: Estimated cumulative density up to each xx point
	"""

	# Parse arguments
	c,xx,xx0,ab = plot_cheb_argparse(1001,*varargin)

	N = len(c)
	txx = np.arccos(xx)
	yy = c[0]*(txx-np.pi)/2
	for idx in np.arange(1,N):
		yy += c[idx]*np.sin(idx*txx)/idx

	yy *= -2/np.pi

	# Plot by default
	if pflag:
		plt.plot(xx0,yy)
		plt.ion()
		plt.show()
		plt.pause(1)
		plt.clf()

	return yy

def plot_chebp(varargin,pflag=True):
	"""
	Given a set of first-kind Chebyshev moments, compute the associated 
	polynomial (*NOT* a density). Output a plot of the polynomial by default.

	Args:
		c: Vector of Chebyshev moments (on [-1,1])
		xx: Evaluation points (defaults to mesh of 1001 pts)
		ab: Mapping parameters (default to identity)
		pflag: Option to output the plot

	Output:
		yy: Polynomial evaluated at xx mesh
	"""

	# Parse arguments
	c,xx,xx0,ab = plot_cheb_argparse(1001,*varargin)

	# Run the recurrence
	kind = 1
	N = len(c)
	P0 = xx*0+1
	P1 = kind*xx
	yy = c[0]+c[1]*xx
	for idx in np.arange(2,N):
		Pn = 2*(xx*P1)-P0
		yy += c[idx]*Pn
		P0 = P1
		P1 = Pn

	# Plot by default
	if pflag:
		plt.plot(xx0,yy)
		plt.ion()
		plt.show()
		plt.pause(1)
		plt.clf()

	return yy

if __name__ == '__main__':
	pass
