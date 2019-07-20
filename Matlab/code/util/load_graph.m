function [A varargout]=load_graph(graphname)
% LOAD_GRAPH Loads a graph from the data directory
%
% load_graph is a helper function to load a graph provided with the
% regardless of the current working directory.  
%
% Example:
%   A = load_graph('cond-mat-2005-fix-cc');




% David F. Gleich
% Copyright, Purdue University, 2013

% History
% 2011-10-02: Initial coding based on load_gaimc_graph
% 2013-01-28: Added option to load coordinates
% 2013-07-21: Removed option for coordinates for spechist
% 2013-07-21: Added option to load eigenvalues

datadir = 'data/rodger'; % The relative directory from NetworkDOS/Matlab/ folder
path=fileparts(mfilename('fullpath'));
path=strrep(path,'code/util',datadir);
smatfile=fullfile(path,[graphname '.smat']);
smatgzfile=fullfile(path,[graphname '.smat.gz']);
if exist(smatfile,'file')
    A=readSMAT(smatfile);
elseif exist(smatgzfile,'file')
    A=readSMAT(smatgzfile);
else
    error('spechist:load_graph','graph %s does not exist',graphname);
end

if nargout > 1
    varargout{1} = load(fullfile(path,[graphname '.smat.normalized.eigs']));
end