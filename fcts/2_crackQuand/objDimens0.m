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
    