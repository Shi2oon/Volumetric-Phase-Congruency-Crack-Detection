function [i_nBC] = inBC(VxLst,sizVOL)
%% ========================================================================
%%       label the voxels didnot existe at any of the previous steps
%% ========================================================================
%% 29/11/2016
%% Yang Chen
%
%       [i_nBC] = inBC(VxLst,sizVOL)
%       ----------------------------------
%

nsteps = length(VxLst);

i_nBC = cell(nsteps,1);

for istep = 1:nsteps
    tic;
    VxLstI = VxLst{istep};
    
    % Did the detected cracks exist before ?  (BC)
    i_BC = false(length(VxLstI),1);
    for ibefore = 1:istep-1
        [VxLstI_BC] = SelectByMask(VxLstI,VxLst{ibefore},sizVOL);
        i0 = ismember(VxLstI,VxLstI_BC);
        i_BC = i_BC | i0;
    end
        clear VxLstI_BC;

    i_nBC{istep} = ~i_BC;
    toc
end
