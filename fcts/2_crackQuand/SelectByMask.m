function [VxSelected] = SelectByMask(VxLstIN,VxLstMask,sizVOL)

cnx = 26;

if iscell(VxLstIN);    VxLstIN = cell2mat(VxLstIN);    end
if iscell(VxLstMask);    VxLstMask = cell2mat(VxLstMask);    end

% construct the binary volume according to the input VxLstIN
VOL0 = false(sizVOL);
VOL0(VxLstIN) = 1;

% connected-component labeling 
CC0 = bwconncomp(VOL0,cnx);
    clear VOL0

% build a label-number volume of the VxLstIN
label3D_0 = labelmatrix(CC0);
    clear CC0

% find the label numbers connecting to the mask
LLst_0 = unique(label3D_0(VxLstMask));
LLst_0(LLst_0==0) = [];

% select the VxLst connecting to the mask
i0 = ismember(label3D_0,LLst_0);
VxSelected = find(i0);

