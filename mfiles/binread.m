%==========================================================================
% BINREAD by André Kurzmann
% read 2D arrays from binary files
%==========================================================================
function A = binread(fname,ny,nx,fmt)

if nargin < 4
    fmt = 'float32';
end

f = fopen(fname,'r');
A = fread(f,inf, fmt);
fclose(f);

if nargin == 2
    nx = length(A)/ny;
    A = reshape(A,ny,nx);
end

if nargin == 3
    A = reshape(A,ny,nx);
end
