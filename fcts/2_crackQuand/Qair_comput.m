function [Qair] = Qair_comput(varargin)
%%%% Compute the air propotion for crack voxels              %%%%
%%%%  - variable gray level gap for different radial postion %%%%
%%%%  - constant gray level gap for all crack voxels         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [Qair] = Qair_comput(Vres,Vgap,p0)
%   ----------------------------------
%  Inputs:
%       > Vres: the gray-levels of crack voxels in the subtracted image
%       > Vgap: the contrast gap between solid and pore
%       > p0:   the maximum interger value of considered tomographic volume
%       > method: the method used for image subtraction
%
%  Outputs:
%       > Qair: the air fraction within crack voxels
%
% 06/04/2016 Y.Chen
% 15/11/2016 Y.Chen
%

narginchk(3,4);

Vres   = varargin{1};
Vgap   = varargin{2};
p0     = varargin{3};
if nargin==4
    method = varargin{4};
else
    method = 'half';
end

if strcmp(method,'half')
    Qair = (single(p0)-single(Vres).*2)./single(Vgap);
    fprintf('formula used: (g-f+255)/2');
elseif strcmp(method,'entire')
    Qair = (single(ceil(p0/2))-single(Vres))./single(Vgap);
    fprintf('formula used: g-f+128');
else
    error('input method invalid');
end

end

