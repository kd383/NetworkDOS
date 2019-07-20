function [t,d,w] = tricnt(G,d,matlabbgl)
%TRICNT Count number of triangles per vertex in a simple, undirected graph
%
%   T = TRICNT(G) takes a sparse adjacency matrix G and computes the number
%   of triangles per vertex. Note taht there is no error checking on G. It
%   is up to the user to ensure that G is symmetric, has only 0/1 entries
%   (but *not* binary), and has no entries on the diagonal.
%
%   T = TRICNT(G,D) takes a second argument which is the degree per vertex
%   and does not recalculate it.
%
%   T = TRICNT(G,D,true) uses the clustering_coefficients from MATLAB_BGL.
%   This assumes that this package is installed and in the path.
%
%   [T,D,W] = TRICNT(G) also returns the degree and number of wedges per
%   vertex.
%
%   NOTE: This is an interface to the MEX function provided by
%   tricnt_mex.c, unless the clustering_coefficients function from
%   MATLAB_BGL is used.
% 
%   See also CCPERDEG.
%
% Tamara G. Kolda, Ali Pinar, and others, FEASTPACK v1.1, Sandia National
% Laboratories, SAND2013-4136W, http://www.sandia.gov/~tgkolda/feastpack/,
% January 2014  

%% License
% Copyright (c) 2014, Sandia National Laboratories
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:  
%
% # Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer. 
% # Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.  
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.          
%
%
% Sandia National Laboratories is a multi-program laboratory managed and
% operated by Sandia Corporation, a wholly owned subsidiary of Lockheed
% Martin Corporation, for the U.S. Department of Energy's National Nuclear
% Security Administration under contract DE-AC04-94AL85000. 


if ~exist('matlabbgl','var')
    matlabbgl = false;
end

if ~exist('d','var') || isempty(d)
    d = full(sum(G,2));
end

w = d.*(d-1)/2;

if (matlabbgl)
    
    if ~exist('clustering_coefficients.m','file')
        error('Must install MATLAB_BGL toolbox');
    end
    options.undirected = 1;
    options.unweighted = 1;
    cc = clustering_coefficients(G,options);
    t = round(w.*cc);
   
else
    
    t = tricnt_mex(G);

end

