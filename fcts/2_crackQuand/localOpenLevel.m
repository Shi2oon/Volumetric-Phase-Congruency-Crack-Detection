function [OpenLv] = localOpenLevel(Vcrk,Vgap,p0,method)
%% compute the local opening level within each crack voxel
%% =======================================================
% 
%       [OpenLv] = localOpenLevel(VxLst,Vsolid,Vgap,method);
%
%
% 01/12/2016  Yang Chen
%


OpenLv = Qair_comput(Vcrk,Vgap,p0,method);

OpenLv(OpenLv>1)=1; OpenLv(OpenLv<0)=0;


