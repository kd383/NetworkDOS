#!/usr/bin/env python
"""
This module is a function to load graph matrix from Matlab file

8/1/2016 - Created
8/16/2016 - Commented
11/25/2016 - Tested
"""

import scipy.io as sio
import numpy as np
import pandas as pd

def load_graph(graphname,path='./data/',mname='A'):
	"""
	Load a graph matrix from Matlab file

	Args:
		graphname: name of the .mat file
		path: the folder .mat file is stored in
		mname: name of the Matlab matrix

	Output:
		data[mname]: the graph matrix
	"""

	data=sio.loadmat(path+graphname)
	return data[mname]

def load_graph_txt(graphname,path='./data/'):
    return pd.read_csv(path+graphname,delimiter='\t',lineterminator='\r',dtype=np.int32).values
    #return np.loadtxt(path+graphname,dtype=np.int)
