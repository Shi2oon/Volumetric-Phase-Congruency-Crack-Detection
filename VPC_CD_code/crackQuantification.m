% clear all;close all

addpath fcts
addpath fcts/2_crackQuand
%mex('fcts\eig3volume_YC.c')

%% data reading ***********************************************************
% read the binary image (cracks labelled as 1)
fname = '../Rdy_for_Avizo/cast_iron/450N_190k_rdy.raw';
% ny = 1601;
% nx = 1645;
% nz = 400;
ny = 1600;
nx = 1600;
nz = 400;

Vbw = lect_raw(fname,nx,ny,nz,5);
% SE = strel('disk',10);tic
% Vbw = imerode(imdilate(Vbw,SE),SE);toc

% read the grey-level image (subtracted image, output of CMV3D)
% fname = '../clean_results/bigvolumeIS.raw';
fname = '../Results_CMV3D/cast_iron/450N/190k/190k_IS.raw';
ny = 1600;
nx = 1600;
nz = 400;
V = lect_raw(fname,nx,ny,nz,5);
% V = 255 - V;
%     %clean the deformed image (remove the grey levels other than cracks)
%     V(~Vbw) = 0;
    

% fname = '../clean_results/bigvolumeIS_tr.raw';
% ny = 1601;
% nx = 1645;
% nz = 400;
% transformed = lect_raw(fname,nx,ny,nz,5);

% %% crack orientation ******************************************************
% %parameters to be set up
% ww = [9, 9, 9]; %window size
% [nR,nC,nB] = size(V);
% % ==== local inertia tensor for all "white" voxel  => local orientation
% nbloc=ceil(nz/50);
% inc=ceil(nz/nbloc);
% Ixx_tot=cell(nbloc,1);  Iyy_tot=cell(nbloc,1);  Izz_tot=cell(nbloc,1);
% Iyz_tot=cell(nbloc,1);  Ixz_tot=cell(nbloc,1);  Ixy_tot=cell(nbloc,1);
% V1_tot=cell(nbloc,1);  V2_tot=cell(nbloc,1);  V3_tot=cell(nbloc,1);tic
% disp('computing the local vector ')
% tic
% for ibloc=1:nbloc
%     i1 = 1+inc*(ibloc-1);
%     i2 = i1+inc-1;
%     if i2>nz; i2=nz; end
%     disp(num2str([i1,i2,i2-i1+1]));
%     % -----------------------------
%         % Inertia Matrix to the center of box
%         VxLst = find(Vbw(:,:,i1:i2));
%         VxLst = VxLst + nR*nC*(i1-1);
% 
%         [Ixx,Iyy,Izz,Iyz,Ixz,Ixy] = I_local(Vbw,VxLst,ww);
% 
%         Ixx_tot{ibloc}=Ixx;
%         Iyy_tot{ibloc}=Iyy;
%         Izz_tot{ibloc}=Izz;
%         Iyz_tot{ibloc}=Iyz;
%         Ixz_tot{ibloc}=Ixz;
%         Ixy_tot{ibloc}=Ixy;
%     
%     % -----------------------------
%         % eigenvalues and eigenvectors of inertia matrices
%         [Lambda1,Lambda2,Lambda3,Vx,Vy,Vz]=eig3volume_YC(Ixx,Ixy,Ixz,Iyy,Iyz,Izz);
%     
%     V1_tot{ibloc}=Vx;
%     V2_tot{ibloc}=Vy;
%     V3_tot{ibloc}=Vz;
% end
% 
% 
% Ixx = cell2mat(Ixx_tot(:));
% Iyy = cell2mat(Iyy_tot(:));
% Izz = cell2mat(Izz_tot(:));
% Iyz = cell2mat(Iyz_tot(:));
% Ixz = cell2mat(Ixz_tot(:));
% Ixy = cell2mat(Ixy_tot(:));
% Vxyz = [cell2mat(V1_tot(:)), cell2mat(V2_tot(:)), cell2mat(V3_tot(:))];
%  
% toc
% 
% % ==== Label the cracks with their sorting number of volsize
% disp('labeling the connected cracks')
% tic
% CC = bwconncomp(Vbw,26);
% vorig=struct2array(regionprops(CC,'Area'));
% [vsort,isort]=sort(vorig,'descend');
% CC_sort = CC;
% CC_sort.PixelIdxList = CC.PixelIdxList(isort);
% label3d_sort=labelmatrix(CC_sort); %~14%x94Go
% Lsort=label3d_sort(Vbw);
% toc
% 
% % ==== preparation for saving
% VxLst    = find(Vbw);
% [y,x,z]  = ind2sub(size(Vbw),VxLst);
% angX_abs = acos(abs(Vxyz(:,1)));
% angY_abs = acos(abs(Vxyz(:,2)));
% angZ_abs = acos(abs(Vxyz(:,3)));
% 
% 
% % ==== visualization
% numCRK = input('choose a crack using their sorting number [1] ');
% if isempty(numCRK);  numCRK=1;  end
% i0 = Lsort == numCRK;
% [figy,figx,figz] = ind2sub(size(Vbw),VxLst(i0));
% c = angX_abs(i0);
% figure;
% scatter3(figx,figy,figz,ones(nnz(i0),1),c);
% axis equal; caxis([0, pi/2]); title(['crack size: ',num2str(nnz(i0))])
% % 
% % ==== save the result
% fnameS = ['crkOrientation.vtk'];
% saveVolvtk_UnStructGrid(fnameS,uint16([x,y,z]),...
%                                'angR_abs',angX_abs,...
%                                'angT_abs',angY_abs,...
%                                'angZ_abs',angZ_abs,...
%                                'Lsort',Lsort,...
%                                'X',uint16(x),...
%                                'Y',uint16(y),...
%                                'Z',uint16(z));
% % % ==== read the saved result
% % [ptPosI,fieldsI] = readVolvtk_UnStructGrid(fnameS,4);
% % VxLstI = sub2ind([nR,nC,nB],ptPosI(:,2),ptPosI(:,1),ptPosI(:,3));
    
%% quantification: thickness and surface area of the cracks ***************
%clean the data (remove the fake cracks)
Vbw = bwareaopen(Vbw,50);
VxLst = find(Vbw);

%prepare for the quantification procedure
Vgap = 100;
[Qair] = Qair_comput(V,Vgap,255,'half');
% Qair(Qair>1) = 1;
% Qair(Qair<0) = 0;
    %figure;hist(Qair(VxLst),100)
    %check----------------------------------
    isl = 223;
    [figy,figx] = find(Vbw(:,:,isl));
    figure;imshow(Qair(:,:,isl))
    hold on;plot(figx,figy,'.')
    figure;imshow(V(:,:,isl))
    hold on;plot(figx,figy,'.')
    %---------------------------------------
%quantification procedure
ww = [15,15,15]; %window size
L = 12; %the range of Qair profile [-L,L]
offset = -min(Qair(:));
% offset = 0;
nblocs = 50;
sizeQair = size(Qair);
QairZ = sizeQair(3);
inc = ceil(QairZ/nblocs);
Wav_tot = [];
WavS_tot = [];
S_tot = [];
dS_tot = [];
Vav_tot = [];
Gav_tot = [];
Nav_tot = [];
over = 0;
tic
for ibloc=1:nblocs
    disp(['--> bloc ',num2str(ibloc),' / ',num2str(nblocs)]);
    i0 = 1+(ibloc-1)*inc;
    i1 = i0+inc-1;
    if i1>QairZ
        i1 = QairZ;
        inc = i1-i0+1;
        over = 1;
    end
    [idX, idY, idZ] = ind2sub(size(Qair), VxLst);
    id = idZ>=i0 & idZ<i1;
    idX = idX(id);
    idY = idY(id);
    idZ = idZ(id);
    tempVxLst = sub2ind(size(Qair), idX, idY, idZ);
    [Wav1,WavS1,S1,dS1,Vav1,Gav1,Nav1] = surfThickness(tempVxLst,Qair,ww,L,1,offset);
    Wav_tot = cat(1, Wav_tot, Wav1);
    WavS_tot = cat(1, WavS_tot, WavS1);
    S_tot = cat(1, S_tot, S1);
    dS_tot = cat(1, dS_tot, dS1);
    Vav_tot = cat(1, Vav_tot, Vav1);
    Gav_tot = cat(1, Gav_tot, Gav1);    %le 1 a changer
    Nav_tot = cat(1, Nav_tot, Nav1);    %idem
    if over == 1
        break
    end
end
toc

%verification by visualization
incv=100;
figure;plot3(Gav_tot(1:incv:end,1),Gav_tot(1:incv:end,2),Gav_tot(1:incv:end,3),'or')
hold on;quiver3(Gav_tot(1:incv:end,1),Gav_tot(1:incv:end,2),Gav_tot(1:incv:end,3),...
                Nav_tot(1:incv:end,1),Nav_tot(1:incv:end,2),Nav_tot(1:incv:end,3),5)
axis equal;xlabel('x');ylabel('y');zlabel('z')




%     %save
%     fnameS = 'dS.bin'; saveBin(dS_tot,fnameS);
%     fnameS = 'Vav.bin';saveBin(Vav_tot,fnameS);
%     fnameS = 'Gav.bin';saveBin(Gav_tot,fnameS);
%     fnameS = 'Nav.bin';saveBin(Nav_tot,fnameS);
%     %read
%     dS = readBin('dS.bin');
%     Vav = readBin('Vav.bin');
%     Gav = readBin('Gav.bin');
%     Nav = readBin('Nav.bin');
    


% lstX = unique(round(Gav1(:,1)));
% lstY = unique(round(Gav1(:,2)));
% lstZ = unique(round(Gav1(:,3)));
% [y,x,z] = ndgrid(lstY,lstX,lstZ);


% Qair(Qair<0) = 0;
% Qair(Qair>1) = 1;
% [Wav2,WavS2,S2,dS2,Vav2,Gav2,Nav2] = surfThickness_2(VxLst,Qair,ww,L,1);

%verification by visualization
% inc=10;
% figure;plot3(Gav1(1:inc:end,1),Gav1(1:inc:end,2),Gav1(1:inc:end,3),'or')
% hold on;quiver3(Gav1(1:inc:end,1),Gav1(1:inc:end,2),Gav1(1:inc:end,3),...
%                 Nav1(1:inc:end,1),Nav1(1:inc:end,2),Nav1(1:inc:end,3),2)
% axis equal
% 
% figure;plot3(Gav2(1:inc:end,1),Gav2(1:inc:end,2),Gav2(1:inc:end,3),'or')
% hold on;quiver3(Gav2(1:inc:end,1),Gav2(1:inc:end,2),Gav2(1:inc:end,3),...
%                 Nav2(1:inc:end,1),Nav2(1:inc:end,2),Nav2(1:inc:end,3),2)
% axis equal

%% COD and crack tip calculation

Coord = [Gav_tot(:,1),Gav_tot(:,2)];
Value = Vav_tot;
Systeme = [Coord, Value];
% Xpixel = 1645/103;
% Ypixel = 1601/100;
Xpixel = nx/200;
Ypixel = ny/200;
finalSysteme = zeros(200,200);
for i = 1:200
    for j = 1:200
        tempSysteme = Systeme(Systeme(:,1) > i*Xpixel & Systeme(:,1) < (i+1)*Xpixel,:);
        tempSysteme = tempSysteme(tempSysteme(:,2) > j*Ypixel & tempSysteme(:,2) < (j+1)*Ypixel,:);
        meanValue = mean(tempSysteme(:,3));
        if isnan(meanValue) || meanValue < 0
            finalSysteme(j,i) = 0;
        else
            finalSysteme(j,i) = meanValue;
        end
    end
end
% 
% Coord_round = round(Coord);
% id = sub2ind(siz([1,2]),Coord_round(:,2),Coord_round(:,1));
% finalSysteme = zeros(siz([1,2]));
% finalSysteme(id) = Value;
% figure;imagesc(finalSysteme');



% finalSysteme = finalSysteme*2.5e-3;
finalSysteme = finalSysteme*2.5e-3;
filtSigma = 4;
imageFilter = strel('disk',4,0);
finalSystemeNaN = finalSysteme;
finalSystemeNaN(finalSystemeNaN==0) = NaN;
finalSystemef = nanconvn(finalSystemeNaN,imageFilter.Neighborhood);
finalSystemef(isnan(finalSystemeNaN)) = 0;

% finalSystemebw = finalSysteme>0;
% r = 16;
% SE = strel('disk',r,0);
% finalSystemef = imdilate(finalSysteme,SE);
% finalSystemef = imerode(finalSystemef,SE);

figure;surf(finalSystemef); colorbar;
xlabel('y (pixel)'); ylabel('x (pixel)'); zlabel('COD measurement (mm)');
title('COD measurement');
xticks([0 50 100 150 200])
xticklabels([0*8 50*8 100*8 150*8 200*8])
yticks([0 50 100 150 200])
yticklabels([0*8 50*8 100*8 150*8 200*8])
% ylim([0 103]); xlim([0 60]);

sizefSf = size(finalSystemef);
contourS = zeros(1,sizefSf(1));
for k = 1:sizefSf(2)
    P = find(finalSystemef(:,k),1);
    if isscalar(P)
        contourS(1,k) = P;
    elseif isempty(P)
        contourS(1,k) = sizefSf(1);
    else
        contourS(1,k) = 0;
    end
end

contourS = contourS(contourS~=0);
figure;plot(contourS(1,:))
hold on;plot(cracktip(1,:),'r')
xlabel('y (pixel)'); ylabel('x (pixel)');
title('Crack tip position');
xlim([0 200]);ylim([0 200]);
xticks([0 50 100 150 200])
xticklabels([0*8 50*8 100*8 150*8 200*8])
yticks([0 50 100 150 200])
yticklabels([0*8 50*8 100*8 150*8 200*8])
legend('Grey-level method','Phase congruency method');


% id = finalSystemef~=0;
% [figy,figx] = find(id);
% incDVCx = 1645/103;
% incDVCy = 1601/100;
% datamap = [1701 - figy.*incDVCy,1585 - figx.*incDVCx];
% notone = ones(length(figx),1)*8;
% Ly = size(transformed,1);
% figure;
% for i =180:220
%     imshow(transformed(:,:,i))
%     hold on;scatter(1701 - figy.*incDVCx,1585 - figx.*incDVCy,notone,finalSystemef(id));
%     waitforbuttonpress
% end
