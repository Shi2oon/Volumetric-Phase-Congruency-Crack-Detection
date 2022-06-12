function [Lsort,CC_sort,vsort,CC] = detectLsort(VxLstI,sizVOL,CONN)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Lable the connected components in descending order of obj volume
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [Lsort,CC_sort,vsort] = detectLsort(VxLstI,sizVOL,CONN)
%   -------------------------------------------------------
%
%   > Inputs:
%     -------
%       > VxLstI: [vector Nvoxels x1] the list of white voxels
%       > sizVOL: [vector 3x1] the dimension size of the considered volume
%       >   CONN: [scalar] the type of 3D connected neighborhood
%
%   > Outputs:
%     -------
%       >   Lsort: [vector Nvoxels x1] the list of sorted label number of
%                                      the input voxel list
%       > CC_sort: [struct] the sorted connectivity information
%       >   vsort: [vector Nvoxels x1] the list of the amount of connected
%                                      voxels
%
%   written by Yang Chen, April 2016
%

% create a binary volume from the input voxel list
VOL = false(sizVOL);
VOL(VxLstI)=1;

% label the connected components of the binary volume
CC = bwconncomp(VOL,CONN);

% sort the labeled components in descending order of volume
vorig=struct2array(regionprops(CC,'Area'));
[vsort,isort]=sort(vorig,'descend');
CC_sort = CC;
CC_sort.PixelIdxList = CC.PixelIdxList(isort);

% get the sorted label numbers
label3d_sort=labelmatrix(CC_sort); %~14%x94Go
Lsort=label3d_sort(VxLstI);

% optimize the data type of vsort
vsort_max = max(vsort);
if vsort_max <= intmax('uint8')
    vsort = uint8(vsort);
elseif vsort_max <= intmax('uint16')
    vsort = uint16(vsort);
elseif vsort_max <= intmax('uint32')
    vsort = uint32(vsort);
end
