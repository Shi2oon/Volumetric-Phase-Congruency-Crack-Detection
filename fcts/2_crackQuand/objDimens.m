function [dr,dt,dz,center,v] = objDimens(VxLst,siz,Ot,ex,ey,ez)
%% %%%%%%%%%%%%%%%%%% measure the dimension of an object   %%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%    in the tube coordinates (RTZ)     %%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [dt,dt,dz,center,v] = objDimens(VxLst,siz,Ot,ex,ey,ez)
%   ------------------------------------------------------
%
%   Inputs:
%       > VxLst: the voxel list of the considered crack
%       > siz: the volume dimension
%       > Ot,ex,ey,ez: tube coordinate bases
%
%   Outputs:
%       > dt: the circumferential width of the considered object
%       > dz: the axial width of the considered object
%       > center: the mass center of the considered object, described in
%                 the tube coordinate (RTZ)
%	> v: the number of voxels of the considered object
%
% written by Yang CHEN 15/04/2016
%

if iscell(VxLst)
    Ncracks = length(VxLst);
    dt = zeros(Ncracks,1);
    dz = dt; dr=dt; v=dt;
    center = zeros(Ncracks,3);
    for icrack=1:Ncracks
        [dr(icrack),dt(icrack),dz(icrack),center(icrack,:),v(icrack)] = ...
                                objDimens0(VxLst{icrack},siz,Ot,ex,ey,ez);
    end

else
    [dr,dt,dz,center,v] = objDimens0(VxLst,siz,Ot,ex,ey,ez);
end

function [dr0,dt0,dz0,center0,v0] = objDimens0(VxLst0,siz,Ot,ex,ey,ez)
    % from voxel list to image coordinates (XYZ)
    [y,x,z] = ind2sub(siz,VxLst0);

    % from image coords (XYZ)  to tube coordinates (RTZ)
    [robj,tobj,zobj] = XYZ2RTZ(Ot,ex,ey,ez,x,y,z);

    % measure the radial width
    dr0 = max(robj) - min(robj);
    
    % measure the circumferential width
    dt0 = max(tobj) - min(tobj);

    if 2*pi-dt0 < pi/100 %if the crack transverses the 2pi angle
        dt0 = 2*pi + max(tobj(tobj<pi/2)) - min(tobj(tobj>pi/2));
    end %hyp: the object has a circum. width smaller than pi

    % measure the axial width
    dz0 = max(zobj) - min(zobj);

    % center of the considered object in the tube coordinate (RTZ)
    center0 = mean([robj,tobj,zobj],1);
    
    % number of voxels of the considered object
    v0 = length(VxLst0);
    
