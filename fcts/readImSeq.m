function [V] = readImSeq(dir0,fmt,varargin)
% read a sequence of images with format of *.tif
%
% [V] = readImSeq(dir0,fmt)
% [V] = readImSeq(dir0,fmt,i0)
% [V] = readImSeq(dir0,fmt,sl0,sl1)
% ---------------------------------------
%
% Inputs:
% -------
%   >    dir0 : path to the files to read
%   >     fmt : extension of image format ('tif','jpg','png',...)
%   >     i0 : reading starts from the i0-th slice
%   >     i1 : reading stops by the i1-th slice
%
% Outputs:
% --------
%   >       V : Grey level table of the whole volume
%
% Yang CHEN, 2018.01.17
%

narginchk(2,4);

% define the names | number of the images to be read
dir0 = [dir0,'/'];
a = dir([dir0,'*.',fmt]);
nsl = length(a);

switch nargin   
    case 2
       i0 = 1;
       i1 = nsl;
    case 3
       i0 = varargin{1};
       i1 = nsl;
    case 4
       i0 = varargin{1};
       i1 = varargin{2};
end

% read images
disp(['--> Reading the image sequence in',dir0])
nn = i1-i0+1;
V = cell(1,1,nn);
for ii = 1:nn
    isl = ii+i0-1;
    fname = [dir0,a(isl).name];
    V{1,1,ii} = imread(fname);
end

% 
V = cell2mat(V);

disp('--> Image sequence reading completed.')

