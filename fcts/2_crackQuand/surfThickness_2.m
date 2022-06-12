function [Wav,WavS,S,dS,Vav,Gav_glob,Nav] = ...
         surfThickness_2(VxLst,VOL_QairI,config,lambda,nsub)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%         Measure the average thickness of a given 3D surface
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  [Wav,S] = surfThickness(VxLst,VOL_Qair,lambda,config,binCapacity,nsub)
%  ----------------------------------------------------------------------
%
%   Inputs:
%   -------
%       >       VxLst: the voxel list describing one surface
%       >   VOL_QairI: the single (double)-valued grayscale volume
%                      corresponding to the considered surface 
%       >      config: the box dimension to be used for computing the local
%                      orientation
%       >      lambda: the interval in the normal direction of the surface
%       >        nsub: discretization degree in the normal direction
%
%
%   Outputs:
%   --------
%       > Wav: average thickness of the suface
%       >   S: estimated surface area of the surface
%
%   written by Yang CHEN 12/07/2016
%

[nR,nC,nB] = size(VOL_QairI); % volume dimension

%%                          determine the local orientation of crack voxels
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ptPosY,ptPosX,ptPosZ]=ind2sub([nR,nC,nB],VxLst);
% Inertia Matrix to the center of box
[Ixx,Iyy,Izz,Iyz,Ixz,Ixy] = I_local(VOL_QairI,VxLst,config);
% eigenvalues and eigenvectors of inertiamatrices
[Lambda1,Lambda2,Lambda3,Vx,Vy,Vz]=eig3volume_YC(Ixx,Ixy,Ixz,Iyy,Iyz,Izz);
VXYZ=[Vx, Vy, Vz];
VXYZ(isnan(Vx),:)=0; %due to Qair(Qair<0)=0; => VxLst~=find(VOLres)
VXYZ(VXYZ(:,3)<0,:) = -VXYZ(VXYZ(:,3)<0,:); %=> Nav=mean(N(i0,1))
VXYZ(VXYZ(:,3)==0&VXYZ(:,2)<0,:) = -VXYZ(VXYZ(:,3)==0&VXYZ(:,2)<0,:);
VXYZ(VXYZ(:,3)==0&VXYZ(:,2)==0&VXYZ(:,1)<0,:) = ...
                            -VXYZ(VXYZ(:,3)==0&VXYZ(:,2)==0&VXYZ(:,1)<0,:);


%%                                            average thickness measurement
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sizSurf = length(VxLst);  %nb of voxels of the given surface


%                        determine the median position for each crack voxel
% =========================================================================

%           interpolation onto the discretized grid in the normal direction
% -------------------------------------------------------------------------
npos_nor = 2*lambda*nsub+1; %discretization 1/nsub voxel
X_nor=zeros(sizSurf,npos_nor);
Y_nor=zeros(sizSurf,npos_nor);
Z_nor=zeros(sizSurf,npos_nor);
w_nor=zeros(sizSurf,npos_nor);
inor = [-lambda:1/nsub:lambda];
for i=1:npos_nor
   X_nor(:,i) = inor(i).*VXYZ(:,1) + ptPosX;
   Y_nor(:,i) = inor(i).*VXYZ(:,2) + ptPosY;
   Z_nor(:,i) = inor(i).*VXYZ(:,3) + ptPosZ;
   
   w_nor(:,i)=interpn(VOL_QairI,Y_nor(:,i),X_nor(:,i),Z_nor(:,i),'linear',0);
end

%                            sum the gray values in every normal direction
%                                       (integrate in the normal direction)
% -------------------------------------------------------------------------
WCrack_median = sum(w_nor,2)./nsub; % x 1/nsub [voxel]


%                 determine the barycenter of the interpolated Qair profile
% -------------------------------------------------------------------------
inor0 = (w_nor*inor'./nsub)./WCrack_median;

X_median = inor0.*VXYZ(:,1)+ptPosX;
Y_median = inor0.*VXYZ(:,2)+ptPosY;
Z_median = inor0.*VXYZ(:,3)+ptPosZ;

clear X_nor Y_nor Z_nor voteId ptPos* wmax


%            Decompose the set of median points into initial voxelized grid
% =========================================================================
X_median_d = round(X_median);
Y_median_d = round(Y_median);
Z_median_d = round(Z_median);


%                         Fusion the rebundant points into one single point
% -------------------------------------------------------------------------
VxLst_median = sub2ind([nR,nC,nB],Y_median_d,X_median_d,Z_median_d);
VxLst_median_uniq = unique(VxLst_median);
nbins = length(VxLst_median_uniq);
Vav=zeros(nbins,1);
Nav=zeros(nbins,3);
Gav = zeros(nbins,3);

for ibin = 1:nbins
    % identity the points beloning to the current bin "ibin"
    i0 = VxLst_median==VxLst_median_uniq(ibin);
    [Y0,X0,Z0] = ind2sub([nR,nC,nB],VxLst_median_uniq(ibin));
        
    % compute the average Qair value of the identified points
    Vav(ibin) = mean(WCrack_median(i0));
    
    % compute the orientation of the surface
    Nav(ibin,:) = mean(VXYZ(i0,:),1);%all VXYZ must have the same signe in z-axis
    
    % compute the average location in the voxel local coordinates
    Gav_glob(ibin,:) = [mean(X_median(i0)),...
                        mean(Y_median(i0)),...
                        mean(Z_median(i0))];
    Gav(ibin,:) = [Gav_glob(ibin,1) - X0,...
                   Gav_glob(ibin,2) - Y0,...
                   Gav_glob(ibin,3) - Z0];
end

% remove the NaN values
i0 = isnan(Vav) | isnan(Gav(:,1));
Nav(i0,:) = [];
Vav(i0) = [];
Gav(i0) = [];
nbins = length(Vav);

% normalize the orientation vectors
Nav = [Nav(:,1)./sqrt(sum( Nav.*Nav,2)),...
       Nav(:,2)./sqrt(sum(Nav.*Nav,2)),...
       Nav(:,3)./sqrt(sum(Nav.*Nav,2))];

%                             compute the elementary area over each bin box
% -------------------------------------------------------------------------
xV  = ones(nbins,1)*[-0.5 0.5 0.5 -0.5 -0.5 0.5 0.5 -0.5];
yV  = ones(nbins,1)*[-0.5 -0.5 0.5 0.5 -0.5 -0.5 0.5 0.5];
zV  = ones(nbins,1)*[-0.5 -0.5 -0.5 -0.5 0.5 0.5 0.5 0.5];
dS = area_element_2(Nav,Gav,xV,yV,zV);

% ---- integration
S = sum(dS(Vav>0));
WavS = Vav'*dS/S;
Wav = mean(Vav);  % 1/S * [sum(Vi*dSi)] = 1/S * [sum(Vi)*sum(dSi)]
                     %                     = 1/N * sum(Vi)
                     %                       <==== if Vi and Si independant


% 
% %% %%%%%%%  verification by visualization
% 
% %% ===== Qair mapping + median surface
% % ---------    Qair mapping
% xmin = min(ptPosX);
% ymin = min(ptPosY);
% zmin = min(ptPosZ);
% x = ptPosX-xmin+1;
% y = ptPosY-ymin+1;
% z = ptPosZ-zmin+1;
% visVOL = zeros(max(y),max(x),max(z),'single');
% vxlstVIS = sub2ind(size(visVOL),y,x,z);
% visVOL(vxlstVIS) = VOL_QairI(VxLst);
% 
% saveVolvtk(visVOL,'verfifCrackQair.vtk','single','point');
% 
% % ----------    median surface
% xM = round(X_median-xmin+1);
% yM = round(Y_median-ymin+1);
% zM = round(Z_median-zmin+1);
% 
% visVOL = zeros(max(y),max(x),max(z),'single');
% vxlstVIS = sub2ind(size(visVOL),yM,xM,zM);
% [vxlstVIS,ia,ic] = unique(vxlstVIS);
% visVOL(vxlstVIS) =  VOL_QairI(VxLst(ia));
% 
% saveVolvtk(visVOL,'verfifCrackMedian.vtk','single','point');

