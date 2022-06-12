function [varargout] = labellist_User(CC,varargin)
%LABELMATRIX_v Create user defined label matrix from BWCONNCOMP structure.
%   Lv = LABELLst_User(CC,...) creates a user defined label list from the
%   connected components structure CC returned by BWCONNCOMP.  
%
%   Inputs:
%   -------
%       - CC         : the BWCONNCOMP structure
%       - l1, l2, ...: user defined labels
%                      the input user defined labels can be any number
%   Outputs:
%   -------
%       - L1, L2, ...: list of voxels labeled by user defined labels
%                      the output number is the same as the number of input
%                      labels
%
% more information, see help labelmatrix
% 
% 23/11/2016 Yang Chen
%

checkCC(CC, mfilename);

%  'uint8', 'uint16', and 'uint32' are the only integer types labelmatrix
%  can be if we want to maintain memory efficiency when calling label2rgb
%  on the output of labelmatrix.

nl = nargin-1;

for i = 1:nl
    
    li = varargin{i};
    lmax = max(li);

    if ~isa(li,'single') && ~isa(li,'double') && ~isa(li,'float')
        if lmax <= intmax('uint8')
            dataType = 'uint8';
        elseif lmax <= intmax('uint16')
            dataType = 'uint16';
        elseif lmax <= intmax('uint32')
            dataType = 'uint32';
        else
            dataType = 'double';
        end
    else 
        dataType = class(li);
    end
    
    VxLst = cell2mat(transpose(CC.PixelIdxList));
    npts = length(VxLst);
    Llst = zeros(npts,1,dataType);

    i0 = 0;
    for k = 1 : CC.NumObjects
        npts_k = length(CC.PixelIdxList{k});
        i1 = i0+npts_k;
        Llst(i0+1:i1) = li(k);
        i0 = i0+npts_k;
    end
    
    
    [tmp,isort] = sort(VxLst,'ascend');
    Llst = Llst(isort);
    
    varargout{i} = Llst;

end


function checkCC(CC,~)
%CHECKCC validates bwconncomp structure

%   Copyright 2008-2011 The MathWorks, Inc.

if ~isstruct(CC)
    error(message('images:checkCC:expectedStruct'));
end

tf = isfield(CC, {'Connectivity','ImageSize','NumObjects','PixelIdxList'});
if ~all(tf)
    error(message('images:checkCC:missingField'));
end



