function [i_BCorAC] = iBCorAC_new(VxLst,sizVOL)
%% ========================================================================
%%          label the voxels existed at the previous step
%%                                     or will remain at the following step
%% ========================================================================
%% 15/11/2016 Yang Chen
%  02/12/2016 Yang Chen
%
%       [i_BCorAC] = iBCorAC_new(VxLst,sizVOL)
%       ----------------------------------
%

nsteps = length(VxLst);

i_BCorAC = cell(nsteps,1);

for istep = 1:nsteps

    disp(['processing the step ',num2str(istep)])

    VxLstI = VxLst{istep};
    npts = length(VxLstI);
    
    % Did the detected cracks exist before ?  (BC)
    i_BC = false(npts,1);
    
    if istep==1
        i_BC = false(npts,1);
    elseif istep>1
        [VxLstI_BC] = SelectByMask(VxLstI,VxLst{istep-1},sizVOL);
        i_BC = ismember(VxLstI,VxLstI_BC);
    end
    clear VxLstI_BC;

    % Will the detected cracks remain afterwards ?  (AC)
    i_AC = false(npts,1);

    if istep==nsteps
        i_AC = false(npts,1);
    elseif istep<nsteps
        [VxLstI_AC] = SelectByMask(VxLstI,VxLst{istep+1},sizVOL);
        i_AC = ismember(VxLstI,VxLstI_AC);
    end
    clear VxLstI_AC

    %
    i_BCorAC{istep} = i_BC | i_AC;

end
