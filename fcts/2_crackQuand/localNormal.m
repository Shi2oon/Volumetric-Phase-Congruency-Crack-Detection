function [Nxyz] = localNormal(VOLopenLv,VxLst,config)
% Compute the local normal vector of each crack voxel
%
% Nxyz = localNormal(VOLopenLv,VxLst,config)
%
% Inputs:
%	> VOLopenLv: the given binary/gray-level volume in which local inertia tensors are computed
%	> VxLst: the chosen points whose inertia tensors are computed
%	> config: vector of 1x3, defining the neighbor box in which local inertia is computed
%
% Outputs:
%	> Nxyz: Local normal vector of each point in VxLst
%
% Yang CHEN 2016


% inertia tensors
[Ixx,Iyy,Izz,Iyz,Ixz,Ixy] = I_local(VOLopenLv,VxLst,config);

% eigenvalues and eigenvectors of inertia tensors
[Lambda1,Lambda2,Lambda3,Nx,Ny,Nz]=eig3volume_YC(Ixx,Ixy,Ixz,Iyy,Iyz,Izz);


Nxyz=[Nx, Ny, Nz];
Nxyz(isnan(Nx),:)=0; %due to OpenLv(OpenLv<0)=0; => VxLst~=find(VOLres)
Nxyz(Nxyz(:,3)<0,:) = -Nxyz(Nxyz(:,3)<0,:); %=> Nav=mean(N(i0,1))
Nxyz(Nxyz(:,3)==0&Nxyz(:,2)<0,:) = -Nxyz(Nxyz(:,3)==0&Nxyz(:,2)<0,:);
Nxyz(Nxyz(:,3)==0&Nxyz(:,2)==0&Nxyz(:,1)<0,:) = ...
                            -Nxyz(Nxyz(:,3)==0&Nxyz(:,2)==0&Nxyz(:,1)<0,:);