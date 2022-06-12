function [i_nBCandAC] = inBCandAC_new(VxLst,sizVOL)
%% ========================================================================
%%          label the voxels existed at any of the previous steps
%%                        or will remain at any of the afterward steps
%% ========================================================================
%% 29/11/2016 Yang Chen
%  05/12/2016 Yang Chen
%
%       [i_BCorAC] = iBCorAC(VxLst,sizVOL)
%       ----------------------------------
%

nsteps = length(VxLst);

i_nBCandAC = cell(nsteps,1);

for istep = 1:nsteps

    disp(['processing the step ',num2str(istep)])

    VxLstI = VxLst{istep};
    npts = length(VxLstI);
    
    % Did the detected cracks exist before ?  (BC)
    i_BC = false(length(VxLstI),1);

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


    if istep<nsteps
        i_nBCandAC{istep} = ~i_BC & i_AC;
    elseif istep==nsteps
        i_nBCandAC{istep} = ~i_BC;
    end
    
end





