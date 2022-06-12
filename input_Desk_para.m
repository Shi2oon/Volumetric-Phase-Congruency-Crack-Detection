% Input : DVC results
% Output : COD field + Visualization
% This code applies the Phase Congruency (PC) method to the input data in
% order to highlight the cracks (discontinuity). It also applies a snake
% algorithm to the result of PC to segment the crack and get a good cracked
% region representation.
% See explaining document for more details on the different parameters and
% lines the user must change.
% This code was originally written by Ahmet Cinar (Sheffield University).
% and maintained and updated in James Marrow Group (University of Oxford)

% last update 11/03/2020 by Abdalrhaman M. Koko update includes
% 			1. adding a pixel to physical values converter
% 			2. 3D visulation (section 7)
%			3. Meshing (section 11)
%			4. Data saving (section 12)
%			5. Made the PC parameters more autonomous, e.g iteration

close all;      restoredefaultpath;   warning ('off'); 	addpath(genpath(pwd));              
DS = com.mathworks.mde.desk.MLDesktop.getInstance();   	DS.closeGroup('Variables');
clc;clear;      set(0,'defaultAxesFontSize',15);       	set(0,'DefaultLineMarkerSize',12) 
tic;    addpath([pwd '\VPC_CD_code']);
addpath([pwd '\VPC_CD_code\AOSLevelsetSegmentationToolboxM']);
addpath([pwd '\VPC_CD_code\Inpaint_nans']);
addpath([pwd '\VPC_CD_code\textprogressbar']);
addpath([pwd '\VPC_CD_code\MatlabFns\FrequencyFilt']);
addpath([pwd '\VPC_CD_code\MatlabFns\PhaseCongruency']);
addpath([pwd '\VPC_CD_code\HexMeshAbaqus']);
addpath([pwd '\VPC_CD_code\shoemake_3D_v04_07_02']);

%% Input data ==========================================================
Dir.inDir       = pwd;      % file directory
Dir.inFile      = 'B00001.dat';
Dir.Subset      = 96;     	% Subset size   
Dir.Overlap     = 75;         	% Overlap in % 
Dir.PixelSize   = 3.25;     % physical size of pixel        
Dir.Pixel_unit  = 'mm';     % units in 'm' 'mm' '\mum'
Dir.StatUs      = 'physical'; % are data in 'physical' or 'pixel' units
Dir.videos_shot = 0;        % 1 if you want a video, 0 if not    

%% DIC info ============================================================
MainDispl = 'Uz';  % Main displacement direction
SliceDir  = 'Uy';  % Direction of slicing
[eulerAngles,rotCentre,Dir.inDir] = ...
    shoemake_3D_v04_07_02_Abdo(Dir.inDir,Dir.inFile,Dir.Pixel_unit,MainDispl);
Dir.inFile = 'Corr-All.txt';
[indexMD, indexSD] = Index_para(MainDispl,SliceDir);    %index equivalent

%% 1. Load data ===========================================================
fprintf('1st: Loading Data .. ');
[RawDisp,VPC,Ux,Uy,Uz,DU] = ...
    LoadData_para(2,1.6,3,0.02,fullfile(Dir.inDir, Dir.inFile),MainDispl,SliceDir);
Dir.Crack   = [Dir.inDir '\PC Analysis ' MainDispl];	mkdir(Dir.Crack)

switch Dir.StatUs 
	case 'physical'	% if already in physical dim.
		[DU] = inUnit(DU,1,0,1,Dir.Pixel_unit); 
	case 'pixel'		% re-scale units
		[DU] = inUnit(DU,Dir.PixelSize,Dir.Overlap,Dir.Subset,Dir.Pixel_unit); 
end

%% 2. Locate crack ========================================================
fprintf('Done\n2nd: locating the crack .. ');
Ncrackedslice = input('2nd: Choose the number of a cracked slice or 999 if you dont know : ');
[Ncrackedslice] = SelectSlice(VPC,indexSD,Ncrackedslice);  
center = center_para(VPC,indexSD,Ncrackedslice);    close

%% 3. 3D Segmentation =====================================================
fprintf('Done\n3rd: 3D Crack Segmentation .. ');
%Parameters to adapt : ittern, mu
%ittern :   number of loop done to properly find the edges, 
%           make sure it is high enough to detect correctly
%mu :       range allowed between inside and outside crack points, 
%           too low will detect non-existent cracks
%           too high won't allow anything to be detected
ittern = length(VPC);
VPC    = (VPC-min(VPC(:)))/(max(VPC(:))-min(VPC(:)));
propagation_weight = 0.015;              GAC_weight = .02; 
% g = ones(size(V)); % linear diffusion 
g        = ac_gradient_map(VPC,1);       delta_t    = 1; 
mu       = 0.2;                          fileHeader = 'B';
numZeros = 5;                            margin     = 1; 

%% 4. Initial condition ===================================================
fprintf('Done\n4th: Set-up initial conditions .. ');
phi = zeros(size(VPC)); 
phi(center(1)-margin:center(1)+margin,center(2)-margin:center(2)+margin,...
    center(3)-margin:center(3)+margin) = 1;

%% 5. MAIN SNAKE ALGORITHM ================================================
fprintf('Done\n5th: Main Snake algorithm for crack detection .. ');
if Dir.videos_shot == 1     % create and open the video
    out = fullfile(Dir.Crack,[Dir.filelist.name ' movie.avi']); 
    v   = VideoWriter(out,'Uncompressed AVI');      open(v); 
end
                                        
for i = 1:ittern
    phi = ac_hybrid_model(VPC-mu, phi-.5, propagation_weight, ...
        GAC_weight, g, ...
        delta_t, 1);
    if exist('h','var') && all(ishandle(h)), delete(h); end
    iso = isosurface(phi,0);
    h = patch(iso,'edgecolor','r','facecolor','w');  
    axis equal;  view(3); xlabel('x');ylabel('y');zlabel('z')
    set(gcf,'name', sprintf('#iters = %d',i));
    sprint_im=['%0' int2str(numZeros) 'd'];
    fileName1=[fileHeader sprintf(sprint_im,i) '.tiff'];
    fileName1=[fileHeader sprintf(sprint_im,i)];
    drawnow; 
    
    if Dir.videos_shot == 1
        title(num2str(i))
        M(i) = getframe(gcf);
        writeVideo(v,M(i))
    end
end
if Dir.videos_shot == 1
    close(v)        % Close and save the video
end
segSize = size(VPC);
Seg3D   = zeros(segSize);
Seg3D(phi>0)=1;     % Logical map of segmentation

%% 6. Prepare for COD =====================================================
fprintf('Done\n6th: Prepare for COD n .. ');
switch MainDispl
    case 'Ux'
        Ux(isnan(Ux))=0;
        F_OutDisp = Ux;
    case 'Uy'
        Uy(isnan(Uy))=0;
        F_OutDisp = Uy;
    case 'Uz'
        Uz(isnan(Uz))=0;
        F_OutDisp = Uz;
end
F_OutDisp(Seg3D==1)=NaN;

%% 7. re-arrange data and 3D visulation of Dat set ========================
fprintf('Done\n7th: 3D Visualisation (make sure to rotate the graphs) .. ');
DU.Seg3D = Seg3D;               DU.Ncrackedslice = Ncrackedslice;
DU = re_arrangeData(indexSD,DU);

plotAllDis(DU.Phy.X,DU.Phy.Y,DU.Phy.Z,abs(DU.Phy.Ux),abs(DU.Phy.Uy),...
    abs(DU.Phy.Uz),DU.Phy.units,MainDispl)
saveas(gcf,[Dir.Crack '\DVC.fig']);
saveas(gcf,[Dir.Crack '\DVC.tif']); close
% with crack taking out
plotAllDis(DU.Phy.X,DU.Phy.Y,DU.Phy.Z,abs(DU.Phy.Seg_Ux),abs(DU.Phy.Seg_Uy),...
    abs(DU.Phy.Seg_Uz),DU.Phy.units,MainDispl)
saveas(gcf,[Dir.Crack '\PC_DVC.fig']);
saveas(gcf,[Dir.Crack '\PC_DVC.tif']); close

%% 8. COD ================================================================
fprintf('Done\n8th: 2D Visualisation .. ');
lastindex = 6 - (indexMD + indexSD);
LC = size(F_OutDisp,lastindex);
for i = 1:size(F_OutDisp,indexSD)
    switch indexSD
        case 1
            STEP = squeeze(F_OutDisp(i,:,:));
        case 2
            STEP = squeeze(F_OutDisp(:,i,:));
        case 3
            STEP = squeeze(F_OutDisp(:,:,i));
    end
    imagesc(STEP);axis equal;
%     pause;
    try
        [cod{i},coda{i}] = COD_para(STEP,lastindex,indexMD);
        
    catch
         cod{i} = zeros(1,LC);
         coda{i} = zeros(1,LC);
    end
    if length(cod{i})<LC
        cl = cod{i};
        cod{i} = zeros(1,LC);
        cod{i}(1:length(cl))=cl;

        cl = coda{i};
        coda{i} = zeros(1,LC);
        coda{i}(LC-length(cl)+1:end)=cl;
    end
    C(i,1:LC) = cod{i};
    D(i,1:LC) = coda{i};  
end
close; C = C.*DU.Pixel2um;
surf(unique(DU.Phy.Y),unique(DU.Phy.X),C); c=colorbar;
% WhitLine(C) % to add a white line at the crack tip/boundries
xlabel(['X [' DU.Pixel_unit ']']);     ylabel(['Y [' DU.Pixel_unit ']']); 
title('COD measurement');               axis image;   
c.Label.String  = DU.Pixel_unit;
saveas(gcf,[Dir.Crack '\' MainDispl '_COD.fig']);
saveas(gcf,[Dir.Crack '\' MainDispl '_COD.tif']); close

%% 9. post-processing ====================================================
fprintf('Done\n9th: post-processing .. ');
x = unique(DU.Phy.X);   stepX = x(2)-x(1);
y = unique(DU.Phy.Y);   stepY = y(2)-y(1);
z = unique(DU.Phy.Z);   stepZ = z(2)-z(1);

iso = isosurface(phi,0);
pts = iso.vertices;
pts(:,3) = 0;   
plot3(pts(:,1).*stepX,pts(:,2).*stepY,pts(:,3).*stepZ,'*') % crack tip
xlabel(['X [' DU.Pixel_unit ']']);     ylabel(['Y [' DU.Pixel_unit ']']); axis image;   
saveas(gcf,[Dir.Crack '\' MainDispl '_Crack_Tip.fig']);
saveas(gcf,[Dir.Crack '\' MainDispl '_Crack_Tip.tif']);       close

sizPhi = size(phi);
im = zeros(sizPhi([1,2]))+1e3;
idx = sub2ind(sizPhi([1,2]),round(pts(:,2)),round(pts(:,1)));
im(idx) = im(idx)*(-1);
imagesc(unique(DU.Phy.X),unique(DU.Phy.Y),im.*DU.Pixel2um)
xlabel(['X [' DU.Pixel_unit ']']);     ylabel(['Y [' DU.Pixel_unit ']']); axis image;
saveas(gcf,[Dir.Crack '\' MainDispl '_Crack_Shape.fig']);
saveas(gcf,[Dir.Crack '\' MainDispl '_Crack_Shape.tif']);     close
    
IM = cat(3,im,im);
crkTip = isosurface(IM,0);
crkTip = crkTip.vertices;
crkTip = crkTip(1:2:end,:);
DU.Phy.crkTip = crkTip.*stepY;
plot(crkTip(:,2).*stepX,crkTip(:,1).*stepY,'*')
xlabel(['X [' DU.Pixel_unit ']']);     ylabel(['Y [' DU.Pixel_unit ']']); axis image; 
XY = xlim;      xlim([XY(1)-stepX*2 XY(2)+stepX*2]);
YX = ylim;      ylim([YX(1)-stepX*2 YX(2)+stepX*2]);
saveas(gcf,[Dir.Crack '\' MainDispl '_Crack_Outliner.fig']);
saveas(gcf,[Dir.Crack '\' MainDispl '_Crack_Outliner.tif']);  close
    
%% 10. Save and trace the coordinate of the crack tip =====================
fprintf('Done\n10th:Save and trace the coordinate of the crack tip .. ');
sizeC = size(C);
cracktip = zeros(1,sizeC(1));
for k = 1:sizeC(2)
    P = find(C(:,k),1);
    if isscalar(P)
        cracktip(1,k) = P;
    elseif isempty(P)
        cracktip(1,k) = sizeC(1);
    else
        cracktip(1,k) = 0;
    end
end
cracktip = cracktip(cracktip~=0);
fk = ind2sub(size(cracktip),find(cracktip==max(cracktip)));
[~, index] = min(abs(fk-mean(fk)));
DU.Phy.tipCordinate = [cracktip(fk(index)).*stepX fk(index).*stepY];
plot(cracktip(fk(index)).*stepX, fk(index).*stepY,'*b'); hold on
plot(cracktip.*stepX,[1:length(cracktip)].*stepY,'--r'); hold off 
xlabel(['X [' DU.Pixel_unit ']']);     ylabel(['Y [' DU.Pixel_unit ']']); axis image; 
XY = xlim;      xlim([XY(1)-stepX*2 XY(2)+stepX*2]);
YX = ylim;      ylim([YX(1)-stepX*2 YX(2)+stepX*2]);
legend('Crack tip position');
saveas(gcf,[Dir.Crack '\' MainDispl '_Crack_tip.fig']);
saveas(gcf,[Dir.Crack '\' MainDispl '_Crack_tip.tif']);  close

%% 11. Mesh and create abaqus structure ===================================
fprintf('Done\n11th:Hexahedron Meshing and create abaqus structure .. ');
BinaryMask_3D = DU.Seg3D;
% convet zeros to ones
BinaryMask_3D(BinaryMask_3D==1)=2 ;          
BinaryMask_3D(BinaryMask_3D==0)=1;
BinaryMask_3D(BinaryMask_3D==2)=0;          
% create 
HexMeshAbaqus([DU.Phy.X(:)  DU.Phy.Y(:)  DU.Phy.Z(:) DU.Phy.Seg_Ux(:) ...
    DU.Phy.Seg_Uy(:) DU.Phy.Seg_Uz(:)],BinaryMask_3D,Dir.Crack,DU.Pixel_unit);

%% 12. savings ============================================================
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   savefig(FigHandle, fullfile(Dir.Crack, [num2str(iFig), '.fig']));  
% end
fprintf('Done\n12th: Saving .. ');
save([Dir.Crack '\OutPut.mat']); 
save([Dir.Crack '\myDVC.mat'],'DU');
fprintf('Done\n All Completed in %.1f minute(s)',toc/60);