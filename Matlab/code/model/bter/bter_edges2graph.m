function [G,G1,G2] = bter_edges2graph(E1,E2)
%BTER_EDGES2GRAPH Create a graph from edge lists.
%
%   G = BTER_EDGES2GRAPH(E1,E2) returns a sparse adjancency matrix
%   corresponding to the given edge lists produced by BTER. The graph is
%   undirected, unweighted, and has no loops, even if E1 and E2 contain
%   these. 
%
%   [G,G1,G2] = BTER_EDGES2GRAPH(E1,E2) returns the graphs corresponding to
%   Phase 1 and Phase 2 in addition to the combined graph.
%
%   See also BTER, EDGES2GRAPH
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

if isempty(E1)
    nnodes = max(E2(:));
elseif isempty(E2)
    nnodes = max(E1(:));
else
    nnodes = max(max(E1(:)),max(E2(:)));
end

if (nargout < 3)
    G = edges2graph([E1;E2],nnodes);
else    
    G1 = edges2graph(E1,nnodes);    
    G2 = edges2graph(E2,nnodes);    
    G = spones(G1+G2);
    G = spdiags(zeros(nnodes,1),0,G);
end

