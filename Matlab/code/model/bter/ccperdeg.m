function [cd,gcc,info] = ccperdeg(G,varargin)
%CCPERDEG Mean clustering coefficient per degree
%
%   CD = CCPERDEG(G) computes the per-degree clustering coefficient, i.e.,
%   CD(d) is the mean clustering coefficient for nodes of degree d. If bins
%   are used, CD(d) returns the clustering coefficient for the bin
%   containing degree d.
%   
%   [CD,GCC] = CCPERDEG(G) also returns the global clustering coefficient.
%
%   [CD,GCC,INFO] = CCPERDEG(G) also returns additional information.
%
%   [...] = CCPERDEG(G,'param',value accepts parameter-value pairs:
%
%   - 'nsamples'  - Number of samples to use. Set to zero for exact
%                   calcuation. Default: 0
%   - 'bins'      - Specify the degree bins for binned data. Default: []
%   - 'tau'       - Specify tau-value for binning. Default: []
%   - 'omega'     - Specify omega-value for binning. Default: []
%   - 'matlabbgl' - Specify use of MATLAB-BGL clusteringcoefficients
%                   function rather than included code. Default: false
%
%   Note that the 'bins' parameters overrides the 'tau' and 'omega'
%   specifications. Otherwise, both 'tau' and 'omega' must be specified to
%   create bins.
%
%   See also TRICNT, BINDATA.
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

% ** Process inputs
params = inputParser;
params.addParamValue('nsamples', 0);
params.addParamValue('bins',[]);
params.addParamValue('tau', []);
params.addParamValue('omega', []);
params.addParamValue('matlabbgl', false);
params.parse(varargin{:});

nsamples = params.Results.nsamples;
bins = params.Results.bins;
tau = params.Results.tau;
omega = params.Results.omega;
matlabbgl = params.Results.matlabbgl;


% ** Create bins
d = full(sum(G,2));
maxd = max(d);

if isempty(bins)
    if isempty(omega) || isempty(tau)
        bins = (1:(maxd+1))';
    else
        nbins = binlookup(maxd+1,omega,tau); 
        bins = binstart((1:(nbins+1))',omega,tau);
    end
end

% **
if nsamples == 0
    
    [t,d,w] = tricnt(G,d,matlabbgl);             
    [~,binId] = histc(d,bins);
    tf = binId > 0;
    binWedges = accumarray(binId(tf),w(tf));
    nbins = length(binWedges);
    binTriangles = accumarray(binId(tf),t(tf),[nbins 1]);
    cdb = binTriangles ./ max(1,binWedges);   
    gcc = sum(t)/sum(w);

else    
   cdb = ccperdegest(G,bins,nsamples); 
   [~,binId] = histc(d,bins);
   tf = binId > 0;
   w = d.*(d-1)/2;
   binWedges = accumarray(binId(tf),w(tf),size(cdb));
   gcc = (binWedges'*cdb) / sum(binWedges);
   t = [];
   binTriangles = [];
end

[~,binId] = histc(1:maxd,bins);
cd(1:maxd,1) = cdb(binId);

% Shorten the bins array to be the same length as cdb
idx = find(cdb > 0, 1, 'last');
cdb = cdb(1:idx);
bins = bins(1:idx);


% Create info
info.nsamples = nsamples;
info.gcc = gcc;
info.bins = bins;
info.cc_per_bin = cdb;
info.deg_per_vertex = d;
info.wedges_per_vertex = w;
info.tris_per_vertex = t;
info.wedges_per_bin = binWedges;
info.tris_per_bin = binTriangles;
