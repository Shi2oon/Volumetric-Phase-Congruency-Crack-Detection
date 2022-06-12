function [idx_neighb] = ind2indneighb(siz,idx,varargin)
%This function returns the voxel list of the input points neighbors (only
%for 3D application, 2D not yet availabel)  
%
%   [idx_neighb] = ind2indneighb(siz,idx)
%   [idx_neighb] = ind2indneighb(siz,idx,shape,config)
%
%   Inputs:
%       siz: the size of image reference (ex. siz=size(VOL0))
%       idx: the indices of input points
%       varargin: for defining a shape of box to "recognize" the neighbors
%           - shape,config
%               shape: string defining the box shape {'cube'(defaut, with a
%                      dimension of 11x11x11),'ellipsoid','cylinder'} 
%               config: vector defining the shape parameters ([row,col,band])
%                       Attention: the dimensions for the cubic box are
%                       2*config+1 !!!
%
%   Outputs:
%       idx_neighb: a matrix of mxn containing the voxel list of neighbors
%                   for all input points. m=length(idx), n=nb_vox_neighbor
%
% 28/09/2015 (Y.Chen)

[shape,config] = ParseInputs(varargin{:});
% 
%disp(['neighborhood recognizing: ',shape,' ',num2str(2.*config+1)])

% type = ['uint',num2str(ceil(log2(prod(siz))))];
type = ceil(log2(prod(siz)));
if type<=8;
    type = 'uint8';
elseif type>8 && type<=16
    type = 'uint16';
elseif type>16 && type<=32
    type = 'uint32';
elseif type>32 && type<64
    type = 'uint64';
else
    disp 'the neighboring table is too heavy !!'
end
    
        [idxR,idxC,idxB] = ind2sub(siz,idx);

switch shape
    case 'cube'
        sizneighb = [config(1)*2+1,config(2)*2+1,config(3)*2+1];
        idx_neighb = zeros(length(idx),prod(sizneighb),type);
        sizcum = cumprod(siz);
        for i=1:prod(sizneighb)
            [ia,ib,ic] = ind2sub(sizneighb,i);
            ia = ia-config(1)-1;
            ib = ib-config(2)-1;
            ic = ic-config(3)-1;
            idx_neighb(:,i) = idx + ic*sizcum(2) + ib*sizcum(1) + ia;
            % out of limit verification
            outoflimitR = idxR+ia<=0 | idxR+ia>siz(1);
            outoflimitC = idxC+ib<=0 | idxC+ib>siz(2);
            outoflimitB = idxB+ic<=0 | idxB+ic>siz(3);
            outoflimit = outoflimitR | outoflimitC | outoflimitB;
            idx_neighb(outoflimit,i) = 0;
        end
    case 'ellipsoid'
        disp '"ELLIPSOID" shape box not yet available'
    case 'cylinder'
        disp '"CYLINDER" shape box not yet available'
end








%%%
%%% ParseInputs
%%%
function [shape,config] = ParseInputs(varargin)

% Check the number of input arguments.
narginchk(0,2);

% Determine the shape and its configuration from the user supplied string
% and values. 
switch nargin
    case 0 %defaut
        shape = 'cube';
        config = [5 5 5];

    case 2
        shape = varargin{1};
        config = varargin{2};
end
shape = validatestring(shape,{'cube','ellipsoid','cylinder'},mfilename,'SHAPE',1);

