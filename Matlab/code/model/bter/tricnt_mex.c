/*  TRICNT_MEX.C: Computes the number of triangles adjacent to each vertex.  

The code uses full enumeration. Each edge is assigned to its lower degree vertex, 
and each vertex checks wedges formed by edges assigned to itself.   

For computational results for this algorithm, see 
C. Seshadhri, A. Pinar, and T.G. Kolda, 
Triadic Measures on Graphs: The Power of Wedge Sampling, 
Proc. SIAM Data Mining, May 2013. 

Tamara G. Kolda, Ali Pinar, and others, FEASTPACK v1.1, Sandia National
Laboratories, SAND2013-4136W, http://www.sandia.gov/~tgkolda/feastpack/,
January 2014  

** License **
Copyright (c) 2014, Sandia National Laboratories
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:  

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer. 

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.  

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.          

Sandia National Laboratories is a multi-program laboratory managed and
operated by Sandia Corporation, a wholly owned subsidiary of Lockheed
Martin Corporation, for the U.S. Department of Energy's National Nuclear
Security Administration under contract DE-AC04-94AL85000.                                         
*/

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <memory.h>

struct graph /* Stored in MATLAB Compressed Sparse Column format */
{
  int V; /* Number of vertices */
  int E; /* Number of edges */  
  mwIndex *ptr; /* ptr[j] = start of column j within ind array */
  mwIndex *ind; /* ind[ptr[j]] = row index for first nonzero in column j */
};

/*  
Checks if u is adjacent to v in G 
Returns 1 if they are adjacent; and 0 otherwise.  
*/
int check_pair(struct graph *G, mwIndex u, mwIndex v)
{
  int i;
  for (i = G->ptr[u]; i < G->ptr[u+1]; i ++)
  {
    if (G->ind[i] == v)
    {
      return(1);
    }
  }
  return(0);
}

/* 
Marks on triangles  formed by vertices on the list "list"  and vertex r
Inputs: G: graph 
list: list of vertices  that are adjacent to r 
n: length of the "list" 
r: vertex r that is the center of wedgesbeing checked
td: an array that stores the number of triangles adjacent to each vertex; 
the array entries are incremented with the new triagles identified     

Output: cnt:  number of triangles found 

*/
int mark_triangles(struct graph *G, mwIndex *list, int n, int r, double *td)
{
  int i, j, cnt, x, y;
  mwIndex *ptr;

  ptr = G->ptr;
  cnt = 0;
  for(i = 0; i < n; i ++)
  {
    x = ptr[list[i]+1] - ptr[list[i]];
    for(j = i+1; j < n; j ++)   /* Check every pair of vertices on the list "list"  to see  if they form a triangle */ 
    {
      if (x < (ptr[list[j]+1] - ptr[list[j]])) /* enables searching via the shorter list */
      {
        y = check_pair(G, list[i], list[j]);
      }
      else 
      {
        y = check_pair(G, list[j], list[i]);  
      }
      if (y)   /*  increment the counters if a triangle is identified */
      {
        cnt ++;
        td[r] ++;
        td[list[j]] ++;
        td[list[i]] ++;
      }
    }
  }
  return(cnt);
}

/* The main function that  counts all the triangles 
Arguments
- G: input graph [unmodified]
- td: array to be filled in with the number of triangles per vertex
Return value
- Total number of triangles (int)

*/   
int tri_enumerate(struct graph *G, double *td)
{
  mwIndex i, j, N, tcnt, t, *d;
  mwIndex *ptr, *ptr2, *ind;

  N = G->V;
  ptr = G->ptr;
  ind = G->ind;

  /* d[i] is the degree of the ith vertex */
  d = (mwIndex*) malloc( sizeof(mwIndex) * N );

  /* ptr2 will be used to shorten adjacency lists, where each edge is assigned the vertex with a smaller degree */
  ptr2 = (mwIndex*) malloc( sizeof(mwIndex) * (N+1) );

  /* Initialize td to zero */
  memset(td, 0, sizeof(double) * N);

  /* make degree list */
  for(i = 0; i < N; i ++)
  {
    d[i] = ptr[i+1] - ptr[i];
    ptr2[i] = ptr[i+1];
  }

  /* Each vertex is assigned to its vertex with a smaller degree 
  edges assigned to a vetex are moved towards the start of each list  such that 
  neighbors  for which the edges assigned to vertex i are listed in ind[ptr[i]] to ind[ptr2[i]-1]
  Note that ind[ptr[i]] to ind[ptr[i]-1] stil stores all neighbors of the ith vertex
  */
  for (i = 0; i < N; i ++) 
  {
    for (j = ptr[i]; j < ptr2[i]; j ++)
    {
      if ((d[i] > d[ind[j]]) || ((d[i] == d[ind[j]]) && (i > ind[j])))
      {
        ptr2[i] --;

        /* swap */
        t = ind[ptr2[i]];
        ind[ptr2[i]] = ind[j]; 
        ind[j] = t; 

        j--;
      }
    }
  }


  /* Check for triangles centered on each vertex with the edges assigned to it */
  tcnt = 0;  
  for (i = 0; i < N; i ++)
  {
    tcnt += mark_triangles(G, ind + ptr[i], ptr2[i] - ptr[i], i, td);
  }

  free(d);
  free(ptr2);
  return(tcnt);
}

/* ----------------------------------------------------------------------------------
This function provides the interface to Matlab 
To call this function, you need to execute  in Matlab the following
>> mex tricnt_mex.c -largeArrayDims

The matlab function sould be called  as 
>> t = tricnt_mex(G)

G is assumed to be a sparse adjacency matrix for a simple graph.
It returns a vector t, such that t[i] is the number of triangles 
adjacent to the ith vertex.
------------------------------------------------------------------------------------ */
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *dtd;
  struct graph G;

  /* Check inputs */
  if ((nrhs != 1) || (!mxIsSparse (prhs[0])) )
  {
    mexErrMsgTxt ("expects sparse matrix");
  }

  /* Read sparse matrix input */
  G.V = mxGetN (prhs [0]);
  G.E = mxGetNzmax (prhs [0]);
  G.ind =  mxGetIr (prhs[0]);
  G.ptr =  mxGetJc (prhs[0]);

  /* Create array for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(G.V, 1, mxREAL);
  dtd = mxGetPr(plhs[0]);

  /* Compute the number of triangles for each vertex */
  tri_enumerate(&G,dtd);
}





