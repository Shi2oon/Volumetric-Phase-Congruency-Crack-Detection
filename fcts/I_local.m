function [Ixx,Iyy,Izz,Iyz,Ixz,Ixy] = I_local(VOLbw,VxLst,config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the local inertia tensor for chosen points in a given gray-level image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [Ixx,Iyy,Izz,Iyz,Ixz,Ixy] = I_local(VOLbw,VxLst,config)
%   --------------------------------------------------------
%
%   Inputs:
%       > VOLbw: the given binary volume in which local inertia tensors are
%                computed
%       > VxLst: the chosen points whose inertia tensors are computed
%       > config: vector of 1x3, defining the neighbor box in which local
%                inertia is computed
%
%   Output:
%       > Ixx,Iyy,Izz,Iyz,Ixz,Ixy
%
%   24/02/2016 Y.Chen
%   
[nR,nC,nB] = size(VOLbw);

% ---- find the neighbors of all points in VxLst
tab_neigh = ind2indneighb([nR,nC,nB],VxLst,'cube',config);
tab_neigh(tab_neigh==0) = find(VOLbw(:,:,1)==0,1,'first'); % => to "remove" the out-limit parts of local box 

% ---- find the phases of neighbors
v_neigh = single(VOLbw(tab_neigh));

% ---- regular grid for neighbor zone (preparation for integration)
[Y0,X0,Z0] = ndgrid(-config(1):config(1),-config(2):config(2),-config(3):config(3));
X0 = single(X0(:));
Y0 = single(Y0(:));
Z0 = single(Z0(:));

% ---- integration to compute the components of inertia tensor
Ixx = v_neigh*(Y0.^2+Z0.^2);
Iyy = v_neigh*(X0.^2+Z0.^2);
Izz = v_neigh*(X0.^2+Y0.^2);
Iyz = v_neigh*(-Y0.*Z0);
Ixz = v_neigh*(-X0.*Z0);
Ixy = v_neigh*(-X0.*Y0);

% % ---- Inertia Matrix to the CoG of box
nneigh = sum(v_neigh,2);
X_G = (v_neigh*X0) ./ nneigh;
Y_G = (v_neigh*Y0) ./ nneigh;
Z_G = (v_neigh*Z0) ./ nneigh;
Ixx = Ixx + nneigh.*(Y_G.^2 + Z_G.^2);
Iyy = Iyy + nneigh.*(X_G.^2 + Z_G.^2);
Izz = Izz + nneigh.*(X_G.^2 + Y_G.^2);
Iyz = Iyz - nneigh.*(Y_G.*Z_G);
Ixz = Ixz - nneigh.*(X_G.*Z_G);
Ixy = Ixy - nneigh.*(X_G.*Y_G);

