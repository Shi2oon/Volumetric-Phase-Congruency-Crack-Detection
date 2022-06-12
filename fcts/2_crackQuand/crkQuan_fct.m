function [Wav,WavS,S,dS,Welem,Gelem,Nelem] = crkQuan_fct(VxLst,OpenLv,Nxyz,sizVOL,lambda,nsub,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%         Measure the average thickness of a given 3D surface
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  [Wav,WavS,S,dS,Welem,Gelem,Nelem] = crkQuan_fct(VxLst,OpenLV,Nxyz,sizVOL,lambda,nsub)
%       ==> usual utilization
%  [Wav,WavS,S,dS,Welem,Gelem,Nelem] = crkQuan_fct(VxLst,OpenLV,Nxyz,sizVOL,lambda,nsub,'visuSaveActivated')
%       ==> unsusual utilization: 
%              save the local openings, orientations, positions into vti file, for the use of verification 
%  ---------------------------------------------------------------------------------------------------------
%
%   Inputs:
%   -------
%       >    VxLst: voxel list describing one surface
%       >   OpenLv: opening level of each crack voxel 
%       >     Nxyz: Orientation(normal) vector of each crack voxel
%       >   sizVOL: dimension of the considered image (ex. 2048,2048,1120)
%       >   lambda: interval in the normal direction of the surface
%       >     nsub: discretization degree in the normal direction
%       > varargin: 'visuSaveActivated', activate the saving of local parameters
%                   into vti file (only for advanced user)
%
%   Outputs:
%   --------
%       >   Wav: average thickness of the surface computed by mean(openings)
%       >  WavS: average thickness of the surface computed by sum(opening*dS)/S
%       >     S: total surface area of the surface
%       >    dS: elementary surface area in each voxel
%       > Welem: elementary opening in each voxel
%       > Gelem: elementary positions in each voxel
%       > Nelem: elementary orientation (normal vector) in each voxel
%
%   written by Yang CHEN 12/07/2016
%   modified by Yang CHEN 25/11/2017 (comments added)
%

[ptPosY,ptPosX,ptPosZ] = ind2sub(sizVOL,VxLst);

VOL_OpenLv = zeros(sizVOL,'single');
VOL_OpenLv(VxLst) = OpenLv;


npts = length(VxLst);  %nb of voxels of the given surface


%                        determine the median position for each crack voxel
% =========================================================================

% 	    interpolation onto the discretized grid in the normal direction
% -------------------------------------------------------------------------
npos_nor = 2*lambda*nsub+1; %discretization 1/nsub voxel
X_nor=zeros(npts,npos_nor);
Y_nor=zeros(npts,npos_nor);
Z_nor=zeros(npts,npos_nor);
w_nor=zeros(npts,npos_nor);
inor = [-lambda:1/nsub:lambda];
for i=1:npos_nor
   X_nor(:,i) = inor(i).*Nxyz(:,1) + ptPosX;
   Y_nor(:,i) = inor(i).*Nxyz(:,2) + ptPosY;
   Z_nor(:,i) = inor(i).*Nxyz(:,3) + ptPosZ;
   
   w_nor(:,i)=interpn(VOL_OpenLv,Y_nor(:,i),X_nor(:,i),Z_nor(:,i),'linear',0);
end

%                            sum the gray values in every normal direction
%                                       (integrate in the normal direction)
% -------------------------------------------------------------------------
WCrack_median = sum(w_nor,2)./nsub; % x 1/nsub [voxel]


%                 determine the barycenter of the interpolated Qair profile
% -------------------------------------------------------------------------
inor0 = (w_nor*inor'./nsub)./WCrack_median;

X_median = inor0.*Nxyz(:,1)+ptPosX;
Y_median = inor0.*Nxyz(:,2)+ptPosY;
Z_median = inor0.*Nxyz(:,3)+ptPosZ;

clear X_nor Y_nor Z_nor ptPos*


%              regroup the set of median points into initial voxelized grid
% =========================================================================
% X_median_d = round(X_median);
% Y_median_d = round(Y_median);
% Z_median_d = round(Z_median);
X_median_d = floor(X_median)+1; % corrected 02/12/2016
Y_median_d = floor(Y_median)+1; % in order to be coordinated exactly 
Z_median_d = floor(Z_median)+1; %             in the image coordinates



%                         Fusion the rebundant points into one single point
% -------------------------------------------------------------------------
VxLst_median = sub2ind(sizVOL,Y_median_d,X_median_d,Z_median_d);
VxLst_median_uniq = unique(VxLst_median);
nbins = length(VxLst_median_uniq);
Welem=zeros(nbins,1);
Nelem=zeros(nbins,3);
Gav = zeros(nbins,3);
Gelem = Gav;


[Y0,X0,Z0] = ind2sub(sizVOL,VxLst_median_uniq);

X0 = X0-0.5;  % corrected 05/12/2016 in order to be coordinated exactly
Y0 = Y0-0.5;  % in the image coordinates
Z0 = Z0-0.5; 

for ibin = 1:nbins
    % identify the points beloning to the current bin "ibin"
    i0 = VxLst_median==VxLst_median_uniq(ibin);

    % compute an average opening value of the identified points
    Welem(ibin) = mean(WCrack_median(i0));
    
    % compute an average normal vector
    Nelem(ibin,:) = mean(Nxyz(i0,:),1);%all Nxyz must have the same signe in z-axis
    
    % compute an average location w.r. to the voxel local coordinates
    Gelem(ibin,:) = [mean(X_median(i0)),...
                     mean(Y_median(i0)),...
                     mean(Z_median(i0))];
               
    Gav(ibin,:) = [Gelem(ibin,1) - X0(ibin),...
                   Gelem(ibin,2) - Y0(ibin),...
                   Gelem(ibin,3) - Z0(ibin)];
end


% remove the NaN values
i0 = isnan(Welem) | isnan(Gav(:,1));
Nelem(i0,:) = [];
Welem(i0) = [];
Gelem(i0,:) = [];
Gav(i0,:) = [];
nbins = length(Welem);

% normalize the normal vectors
Nelem = [Nelem(:,1)./sqrt(sum(Nelem.*Nelem,2)),...
       Nelem(:,2)./sqrt(sum(Nelem.*Nelem,2)),...
       Nelem(:,3)./sqrt(sum(Nelem.*Nelem,2))];

%                             compute the elementary area over each bin box
% -------------------------------------------------------------------------
xV  = ones(nbins,1)*[-0.5 0.5 0.5 -0.5 -0.5 0.5 0.5 -0.5];
yV  = ones(nbins,1)*[-0.5 -0.5 0.5 0.5 -0.5 -0.5 0.5 0.5];
zV  = ones(nbins,1)*[-0.5 -0.5 -0.5 -0.5 0.5 0.5 0.5 0.5];

dS = area_element_2(Nelem,Gav,xV,yV,zV);

% dS(isnan(dS)) = 0;
% if isempty(dS)
%     S = 0;
%     dS = 0;
%     Welem = 0;
%     Gelem = [0 0 0];
%     Wav = 0;
%     WavS = 0;
%     return
% end


%			   compute the total surface area & average opening
% -------------------------------------------------------------------------
S = sum(dS(Welem>0));
if S==0
    WavS=0;  Wav=0;  S=0;
    dS=0;
    Welem=0;
    Gelem=[0 0 0];
    Nelem=[0 0 0];
    return;
end
    

WavS = Welem'*dS/S;    WavS(isnan(WavS)) = 0;
Wav = mean(Welem);  % 1/S * [sum(Vi*dSi)] = 1/S * [sum(Vi)*sum(dSi)]
                     %                     = 1/N * sum(Vi)
                     %                       <==== if Vi and Si independant

                     
                     
% #########################################################################
                     
%% ==== visualization for checking the quantification procedure

if ~isempty(varargin)
    if strcmp( varargin{1},'visuSaveActivated' )
        visuSaveC = input(['are u sure to save the local opening level ',...
                           'into a vti-file ? and to the projection',...
                           ' positions into a vtk-unstructured-file ? (y/n) '],'s');
    end
    
    if strcmp( visuSaveC,'y' )

        %% ==== OpenLv mapping
        % ---------    Qair mapping
        VOL_OpenLV = zeros(sizVOL,'single');
        VOL_OpenLV(VxLst) = OpenLv;
        
        sauv_vti(VOL_OpenLv,1,1,1,...
                 'crkQuanVerif_OpenLv.vti','float','OpenLv');

        %% ==== projection positions
        
        saveVolvtk_UnStructGrid('crkQuanVerif_OpenLv.vtk',single(Gelem-0.5),...
                                'Nx',Nelem(:,1),...
                                'Ny',Nelem(:,2),...
                                'Nz',Nelem(:,3),...
                                'Welem',Welem,...
                                'x',single(Gelem(:,1)-0.5),...
                                'y',single(Gelem(:,2)-0.5),...
                                'z',single(Gelem(:,3)-0.5) );
    end
    
end
