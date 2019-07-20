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
import numpy as np
import scipy as sp
import time
import pandas as pd

def test_moment_comp(L,n,Z=0,x=0):
    if isinstance(Z,int):
        Z = np.sign(nr.randn(n,20))
    if isinstance(x,int):
        x = nr.randn(10,1)

    L = matrix_normalize(A)
    f = open('./result.txt','a',0)
    start = time.time()
    [c, cs] = moments_cheb_dos(L,n,Z,10)
    end = time.time()
    f.write('DOS time: '+str(end-start)+'\n')
    start = time.time()
    [cl, csl] = moments_cheb_ldos(L,n,Z,10)
    end = time.time()
    f.write('LDOS time: '+str(end-start)+'\n')
    f.close()

if __name__ == '__main__':
    f = open('./result.txt','a',0)
    start = time.time()
    #ind = load_graph_txt('dblp.txt','../data/')
    #ind = np.loadtxt('../data/dblp.txt',dtype=np.int)
    ind = pd.read_csv('../data/facebook.txt',delimiter = ' ',
    lineterminator='\n',dtype=np.int64,header=None).values
    end = time.time()
    print(end-start)
    f.write('Load file time: '+str(end-start)+'\n')
    start = time.time()
    imax = np.max(np.max(ind))
    imap = np.full(imax+1, -1, dtype=np.int)
    k = 0
    for j in np.arange(2):
        for i in np.arange(ind.shape[0]):
            if imap[ind[i,j]]==-1:
                imap[ind[i,j]] = k
                k = k+1
            ind[i,j] = imap[ind[i,j]]
    end = time.time()
    f.write('Reorder index time: '+str(end-start)+'\n')
    start = time.time()
    deg = np.zeros(k)
    for i in np.arange(ind.shape[0]):
        deg[ind[i,0]]+=1
        deg[ind[i,1]]+=1
    end = time.time()
    for i in np.arange(k):
        deg[i] = np.sqrt(deg[i])
    f.write('Degree count time: ' + str(end-start)+'\n')
    start = time.time()
    weight = np.zeros(ind.shape[0])
    for i in np.arange(ind.shape[0]):
        weight[i] = 1.0/deg[ind[i,0]]/deg[ind[i,1]]
    end = time.time()
    f.write('Normalizing time: '+str(end-start)+'\n')
    # ind = np.reshape(ind,(-1,1))
    # _,ind = np.unique(ind,return_inverse=True);
    # ind = np.reshape(ind,(-1,2))
    # ind = np.concatenate((ind, np.fliplr(ind)));
    n = 4039
    start = time.time()
    A = sp.sparse.csr_matrix((weight,(ind[:,0],ind[:,1])), shape = (n,n))
    end = time.time()
    f.write('Form matrix time: '+str(end-start)+'\n')
    start = time.time()
    A = A+ A.T
    end = time.time()
    f.write('Symmetrize time: '+str(end-start)+'\n')
    f.write(str(A.shape[0])+'\n')
    f.close()
    test_moment_comp(A,A.shape[0])
