#!/bin/sh
#
# Download network eigenvalue tarball from David Gleich and 
# rename the resulting directory to data.
#
wget https://www.cs.purdue.edu/homes/dgleich/rodger/graph-eigs-v1.tar.gz
tar -xzf graph-eigs-v1.tar.gz
mv graph-eigs-v1 rodger
rm graph-eigs-v1.tar.gz
