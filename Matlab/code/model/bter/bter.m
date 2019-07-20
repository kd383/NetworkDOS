function [e1,e2,info] = bter(nd,cd,varargin)
%BTER Generate edges for BTER model.
%
%   [E1,E2] = BTER(ND,CD) generates a list of edges for the BTER model for
%   the given degree distribution and clustering coefficients per degree,
%   i.e., ND(d) = the number of nodes of degree d and CD(d) = mean
%   clustering coefficient for degree d. The edge list E = [E1;E2] is a
%   list of edges created by the BTER model from phases 1 and 2,
%   respectively. The edge list will have duplicates, which can be remedied
%   in post-processing.    
%
%   [E1,E2,INFO] = BTER(...) returns extra information.
%
%   [...] = BTER(ND,CD,'param',value,...) takes a sequence of
%   parameter-value pairs to specify the operation of the method
%
%   - 'blowup'      - Multiplier for number of degree-1 nodes. Default: 1.
%   - 'rmloops'     - Remove self-links. Default: true.
%   - 'cleanup'     - Clean-up degree-1 nodes per blowup. Default: true.
%   - 'swap'        - Swap edge endpoints so least is 1st. Default: false.
%   - 'verbose'     - Print information. Default: true.
%   - 'rngseed'     - Seed for the random number generator.
%
%   References:
%   * C. Seshadhri, T. G. Kolda and A. Pinar. Community structure and
%     scale-free collections of Erdös-Rényi graphs, Physical Review E
%     85(5):056109, May 2012. (doi:10.1103/PhysRevE.85.056109)
%   * T. G. Kolda, A. Pinar, T. Plantenga and C. Seshadhri. A Scalable
%     Generative Graph Model with Community Structure,  arXiv:1302.6636,
%     March 2013. (http://arxiv.org/abs/1302.6636)
%
%   See also BTER_EDGES2GRAPH, CCPERDEG
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


% ** Check inputs
if nargin < 1
    error('At least one input is required');
end

if nargout < 1
    error('At least one output is required');
end

% ** Process inputs 
params = inputParser;
params.addParamValue('blowup', 1);
params.addParamValue('swap', false);
params.addParamValue('rmloops', true);
params.addParamValue('cleanup', true);
params.addParamValue('verbose', true);
params.addParamValue('rngseed', []);
params.parse(varargin{:});

beta = params.Results.blowup;
if beta < 1
    error('Must have ''blowup'' >= 1.');
end

swap = params.Results.swap;
rmloops = params.Results.rmloops;
cleanup = params.Results.cleanup;
verbose = params.Results.verbose;
rngseed = params.Results.rngseed;

% Set random number generator
if isempty(rngseed)
    rngseed = rng;
end
rng(rngseed);

% Input checking
dmax = length(nd);
if length(cd) ~= dmax
    error('Inputs ND and CD must be the same length');
end

if any(nd < 0)
    error('Degree distribution cannot be negative');
end
if any(round(nd) ~= nd)
    error('Degree distribution must be integral');
end
if any(cd < 0) || any(cd > 1)
    error('Clustering coefficients must be between 0 and 1');
end

% ** Setup 
tic;
[id,wd,ndfill,rdfill,ig,wg,bg,ng] = bter_setup(nd, cd, beta);
time_setup = toc;

if (verbose)
    fprintf('--- BTER HPC Set-up ---\n');
    fprintf('Desired number of nodes: %d\n', sum(nd));
    fprintf('Desired number of edges: %.0f\n', sum(nd.*(1:dmax)')/2);
    fprintf('Multiplier to degree-1 nodes: %g\n', beta);
    fprintf('Maximum degree: %d\n', dmax);
    fprintf('Number of groups: %d\n', length(ig));
    fprintf('Number of blocks: %d\n', sum(bg));   
    fprintf('Phase 1 total weight: %.0f\n', sum(wg));
    fprintf('Phase 2 total weight: %.0f\n', sum(wd));
    fprintf('Time for setup (sec): %.2f\n', time_setup);
end


% ** Determine number of samples per phase 
tic;

w1 = sum(wg);
w2 = sum(wd);
w = w1+w2;
nsmp = round(w);
r = rand(nsmp,1);
s1 = sum(r < w1/w);
s2 = nsmp - s1;

time_split = toc;

if (verbose)
    fprintf('Determined phase for %d edges in %f seconds\n', nsmp, time_split);
end

% ** Phase 1 Samples 
tic;

grp_smp = random_sample(wg, s1);
blk_r = rand(s1,1);
blk_b = bg(grp_smp);
blk_i = ig(grp_smp);
blk_n = ng(grp_smp);
e1_shift = blk_i + floor(blk_r .* blk_b) .* blk_n;
e1_r = rand(s1,2);
e1(:,1) = floor(e1_r(:,1) .* blk_n) + e1_shift;
e1(:,2) = floor(e1_r(:,2) .* (blk_n - 1)) + e1_shift;
e1(:,2) = e1(:,2) + (e1(:,2) >= e1(:,1));

time_p1 = toc;

if (verbose)
    fprintf('Sampled %d edges for phase 1 in %f seconds\n', s1, time_p1);
end

% ** Phase 2 Samples 
tic;

% Setup
idfill = id;
idbulk = id + ndfill;
ndbulk = nd - ndfill;
ndbulk(1) = 0;
% Sample
d_smp = random_sample(wd, 2*s2);
d_smp = reshape(d_smp, s2, 2);
r = rand(s2,2);
tf_fill = r < rdfill(d_smp);
e2_shift_fill = idfill(d_smp);
e2_sz_fill = ndfill(d_smp);
e2_shift_bulk = idbulk(d_smp);
e2_sz_bulk = ndbulk(d_smp);

r = rand(s2,2);
e2_fill = e2_shift_fill + floor(r .* e2_sz_fill);
e2_bulk = e2_shift_bulk + floor(r .* e2_sz_bulk);
e2 = tf_fill .* (e2_fill) + ~tf_fill .* (e2_bulk);

time_p2 = toc;

if (verbose)
    fprintf('Sampled %d edges for phase 2 in %f seconds\n', s2, time_p2);
end

% ** Remove loops (Phase 2 only) 
if rmloops
    tic;
    tf = (e2(:,1) == e2(:,2));
    e2 = e2(~tf,:);
    time_rmloops = toc;   
    if (verbose)
        fprintf('Removed %d loops in %f seconds\n', sum(tf), time_rmloops);
    end
end

% ** Reorder edges so lower index is first 
if swap
    tic;
    e1 = bter_swap(e1);
    e2 = bter_swap(e2);
    time_swap = toc;
    if (verbose)
        fprintf('Swapped low degree first in %f seconds\n', time_swap);
    end
end

% ** Reindexing "degree-1" vertices to only keep those that have a link 
if (beta > 1) && (cleanup)
    tic;
    idx = id(1); % Index of first degree-1 node
    dtrue = accumarray(e2(:),1);
    tf = ones(size(dtrue)) > 0;
    tf(idx:end) = dtrue(idx:end) > 0;
    n_old = length(dtrue);
    n_new = sum(tf);
    old2newidx = zeros(length(dtrue),1);
    old2newidx(tf) = (1:n_new)';
    e2 = old2newidx(e2);
    time_cleanup = toc;
    if (verbose)
        fprintf('Removed %d spurious degree-1 nodes in %f seconds\n', ...
            n_old - n_new, time_cleanup);
    end
end

% ** Finishing up 
if (verbose)
    fprintf('--- BTER HPC Complete ---\n');
end

% Fill in INFO
if nargout >= 3
    S = whos;
    info = struct;
    for i = 1:length(S);
        var = S(i).name;
        if (var(1) == 'e') || strcmp(var(1:min(3,end)),'blk') || ...
                strcmp(var, 'old2newidx') || strcmp(var, 'd_smp') || ...
                strcmp(var, 'grp_smp') || strcmp(var, 'tf') || ...
                strcmp(var, 'dtrue') || strcmp(var, 'r') || ...
                strcmp(var,'tf_fill')
            continue;
        end
        eval(['info.' var ' =  ' var ';']);
    end
end

function edges = bter_swap(edges)
%BTER_SWAP Given a list of edges, make lower index first.

tf = edges(:,1) > edges(:,2);
tmp = edges(tf,2);
edges(tf,2) = edges(tf,1);
edges(tf,1) = tmp;

function [id,wd,ndfill,rdfill,ig,wg,bg,ng] = bter_setup(nd, cd, beta)
%BTER_SETUP Creates data for generating BTER graphs.

% Compute maximum degree
dmax = length(nd);

% Set up arrays (max # groups = dmax)
id = zeros(dmax,1); % i_d 
wd = zeros(dmax,1); % w_d
rdfill = zeros(dmax,1); % r^{\rm fill}_d
ndfill = zeros(dmax,1); % n^{\rm fill}_d
wg = zeros(dmax,1); % w_g
ig = zeros(dmax,1); % i_g
bg = zeros(dmax,1); % b_g
ng = zeros(dmax,1); % n_g

% Index of first node for each degree. 
% Degree 1 vertices are numbered last.
tmp = cumsum(nd(2:end));
id(2) = 1;
id(3:end) = tmp(1:end-1)+1;
id(1) = tmp(end)+1;

% Compute number of nodes with degree greater than d
tmp = flipud(cumsum(flipud(nd)));
ndprime = zeros(dmax, 1);
ndprime(2:end-1) = tmp(3:end);

% Handle degree-1 nodes
ndfill(1) = nd(1) * beta;
wd(1) = 0.5 * nd(1);
rdfill(1) = 1;

% Main loop
g = 0;
nfillblk = 0;
intdeg = 0;
for d = 2:dmax
    
    if nfillblk > 0
        ndfill(d) = min( nfillblk, nd(d) );
        nfillblk = nfillblk - ndfill(d);
        wdfilltmp = 0.5 * ndfill(d) * (d - intdeg);
    else
        ndfill(d) = 0;
        wdfilltmp = 0;
    end
    
    ndbulktmp = nd(d) - ndfill(d);
    
    if ndbulktmp > 0
        g = g + 1;
        ig(g) = id(d) + ndfill(d);
        bg(g) = ceil(ndbulktmp / (d+1));
        ng(g) = d+1;
        if (bg(g) * (d+1)) > (ndprime(d) + ndbulktmp)
            if bg(g) ~= 1
                error('Last group has more than 1 block');
            end
            ng(g) = ndprime(d) + ndbulktmp;
        end
        rho = nthroot(cd(d), 3);
        intdeg = (ng(g) - 1) * rho;
        wdbulktmp = 0.5 * ndbulktmp * (d - intdeg);        
        wg(g) = bg(g) * 0.5 * ng(g) * (ng(g) - 1) * log (1/(1-rho));  
        nfillblk = bg(g) * ng(g) - ndbulktmp;
    else
        wdbulktmp = 0;
    end
    
    wd(d) = wdbulktmp + wdfilltmp;
    if (wd(d) > 0)
        rdfill(d) = wdfilltmp / wd(d);
    else
        rdfill(d) = 0;
    end        
end

% Shorten the group arrays
ig = ig(1:g);
wg = wg(1:g);
bg = bg(1:g);
ng = ng(1:g);



