function A = ImFmtConvert(A,fmt,varargin)
%convert an image to another format ('uint8' or 'uint16')
%
%   A = ImFmtConvert(A,fmt);
%   A = ImFmtConvert(A,fmt,[vmin,vmax]);
%   ------------------------------------
%
%   Inputs:
%       >  A : the image to be converted
%       > fmt: target format (choose between 'uint8' and 'uint16')
%       > [vmin,vmax]: optional, specifying the min and max for format
%                      converting. If absent, vmin/max=min/max(A(:));
%
%   Outputs:
%       >  A : the image of new format
%
% Yang CHEN 2018.04.04
%

% prepare 0 : the ends of interval
if nargin==2
    vmin = double(min(A(:)));
    vmax = double(max(A(:)));
elseif nargin>2
    vmin = varargin{1}(1);
    vmax = varargin{1}(2);
end


% format converting
interv0 = vmax - vmin;
if strcmp(fmt,'uint8')
    interv1 = 256;
    A = (A-vmin) .* (interv1/interv0) ;
    A = uint8(A);
elseif strcmp(fmt,'uint16')
    interv1 = 65536;
    A = (A-vmin) .* (interv1/interv0) ;
    A = uint16(A);
end
