function [VxLstI] = updatCrkVox(VxLstI,i_BCorAC_I,VOLres,method,threVOL0)
%% ========================================================================
%%          update the crack voxels for quantification
%%              - remove the "closed"     (Qair<=0)
%%              - remove the small cracks (<threVOL0)
%% ========================================================================
%% 16/11/2016
%% Yang Chen
%
%       VxLstI = updatCrkVox(VxLstI,i_BCorAC_I,VOLres,method,threVOL0)
%       --------------------------------------------------------------
%

sizVOL = size(VOLres);
p0 = double(intmax(class(VOLres)));

% select the BC-or-AC remaining cracks
VxLstI = VxLstI(i_BCorAC_I);

% remove the "closed"
Vres = VOLres(VxLstI);
if strcmp(method,'half')
    i0 = (p0 - single(Vres).*2) <=0;
elseif strcmp(method,'entire')
    i0 = (p0/2 - single(Vres)) <=0;
end
VxLstI(i0) = [];

% relabel the "opening"
[LsortI,CC_sortI,vsortI] = detectLsort(VxLstI,sizVOL,26);
VxLstI = CC_sortI.PixelIdxList;

% remove the small cracks
i0 = vsortI<threVOL0;
VxLstI(i0) = [];
VxLstI = VxLstI';
        clear LsortI CC_sortI vsortI;