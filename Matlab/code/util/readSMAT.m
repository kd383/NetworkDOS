function A = readSMAT(filename)
% readSMAT reads an indexed sparse matrix representation of a matrix
% and creates a MATLAB sparse matrix.
%
%   A = readSMAT(filename)
%   filename - the name of the SMAT file
%   A - the MATLAB sparse matrix
%

% David Gleich
% Copyright, Stanford University, 2005-2013

% History
% :2013-07-21: Added gzip read ability

if (~exist(filename,'file'))
    error('readSMAT:fileNotFound', 'Unable to read file %s', filename);
end

if (exist(strcat(filename, '.info'), 'file'))
    s = load(filename);
    mdata = load(strcat(filename, '.info'));
    ind_i = s(:,1)+1;
    ind_j = s(:,2)+1;
    val = s(:,3);
    A = sparse(ind_i,ind_j,val, mdata(1), mdata(2));
    return;
end;

[pathstr,name,ext] = fileparts(filename);
if ~isempty(strfind(ext,'.gz'))
    [m n i j v] = readSMATGZ(filename);

    A = sparse(i,j,v,m,n);
    return;
end;


s = load(filename,'-ascii');
m = s(1,1);
n = s(1,2);
try
    ind_i = s(2:length(s),1)+1;
    ind_j = s(2:length(s),2)+1;
    val = s(2:length(s),3);
    clear s;
    A = sparse(ind_i,ind_j,val, m, n);
catch
    fprintf('... trying block read ...\n');
    blocksize = 1000000;
    curpos = 2;
    blocknum = 1;
    nzleft = s(1,3);
    A = sparse(m,n);
    while (nzleft > 0)
        curblock = min(nzleft, blocksize);
        fprintf('block %i (%i - %i)\n', blocknum, curpos-1, curpos+curblock-2);
        curpart = curpos:(curpos+curblock-1);
        ind_i = s(curpart,1)+1;
        ind_j = s(curpart,2)+1;
        val = s(curpart,3);
        A = A + sparse(ind_i, ind_j, val, m, n);
        
        nzleft = nzleft - curblock;
        curpos = curpos + curblock;
        blocknum = blocknum + 1;
    end
    
end


function S = merge_structs(A, B)
% MERGE_STRUCTS Merge two structures.
%
% S = merge_structs(A, B) makes the structure S have all the fields from A
% and B.  Conflicts are resolved by using the value in A.
%

%
% merge_structs.m
% David Gleich
%
% Revision 1.00
% 19 Octoboer 2005
%

S = A;

fn = fieldnames(B);

for ii = 1:length(fn)
    if ~isfield(A, fn{ii})
        S.(fn{ii}) = B.(fn{ii});
    end
end

function [m n is js vs] = readSMATGZ(filename)

% Create Java input stream from the gzipped filename.
fileInStream = [];
try
   fileInStream = java.io.FileInputStream(java.io.File(filename));
catch exception
   % Unable to access the gzipped file.
   if ~isempty(fileInStream)
     fileInStream.close;
   end
   eid = sprintf('MATLAB:%s:javaOpenError',mfilename);
   error(eid,'Could not open file "%s" for reading.',filename);
end

% Create a Java GZIP input stream from the file input stream.
try
   gzipInStream = java.util.zip.GZIPInputStream( fileInStream );
catch exception
   % The file is not in gzip format.
   if ~isempty(fileInStream)
     fileInStream.close;
   end
   eid = sprintf('MATLAB:%s:notGzipFormat',mfilename);
   error(eid,'File "%s" is not in GZIP format.',filename);
end

gzipChannel = java.nio.channels.Channels.newChannel(gzipInStream);
bufsize = 16*1024*1024; % 16MB buffer
B = java.nio.ByteBuffer.allocate(bufsize);

% read a chunk
B.clear();
len = gzipChannel.read(B);
str = char(B.array()');
filedone = false;
if len < length(str)
    str = str(1:len);
    filedone = true;
end
% find the 
if ~filedone
    lastnl = find_last_newline(str);
    strextra = str(lastnl+1:end);
    str = str(1:lastnl);
else 
    strextra = [];
end

% read the header
[s,pos] = textscan(str, '%f',3);
str = str(pos+1:end); % truncate the string to what was used.

m = s{1}(1);
n = s{1}(2);
nnz = s{1}(3);

is = zeros(nnz,1);
js = zeros(nnz,1);
vs = zeros(nnz,1);

curind = 1;

while 1
    
    [buf,pos] = textscan(str,'%f %f %f');
    indi = buf{1};
    indj = buf{2};
    vals = buf{3};
    
    numnewnz = numel(indi);
    is(curind:curind+numnewnz-1) = indi+1;
    js(curind:curind+numnewnz-1) = indj+1;
    vs(curind:curind+numnewnz-1) = vals;
    curind = curind + numnewnz;
    
    str = str(pos+1:end);
    assert(isempty(str)); % we should have used all of the string
    
    if filedone
        % if filedone was set, we are done!
        break
    end
    
    % read the next file chunk
    B.clear();
    len = gzipChannel.read(B);
    str = char(B.array()');
    filedone = false;
    if len < length(str)
        str = str(1:len);
        filedone = true;
    end
    % add in the old extra
    if ~isempty(strextra)
        str = [strextra str];
    end
    if ~filedone
        % find the last newline
        lastnl = find_last_newline(str);
        strextra = str(lastnl+1:end);
        str = str(1:lastnl);
    else 
        strextra = [];
    end
end

if ~isempty(fileInStream)
 fileInStream.close;
end


function ind=find_last_newline(str)

ind = [];
for i=length(str):-1:1
    if str(i) == 10 || str(i) == 13 % these are the newline characters
        ind = i;
        break
    end
end