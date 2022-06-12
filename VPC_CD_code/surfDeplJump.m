function [] = surfDeplJump(Gav,Nav, lambda)
     
fname = '..\DVC_sub\300N_292k-2000N_290k_bigvol_32laststep_bigvol/B00002.dat';

fid=fopen(fname,'r');
ligne = fgetl(fid);
ligne = fgetl(fid);
ligne = fgetl(fid);
fclose(fid);
nx = str2num(ligne(strfind(ligne,'I=')+2:strfind(ligne,', J=')));
ny = str2num(ligne(strfind(ligne,'J=')+2:strfind(ligne,', K=')));
nz = str2num(ligne(strfind(ligne,'K=')+2:strfind(ligne,', F=')));
siz_DVC = [nx ny nz]; %...lire dans *.dat
incDVC = 32*(1-0.75);

A = dlmread(fname,' ',3,0); %read the file
U_pos = A(:,1:3)/(2.5e-03);
U = A(:,4:6);
id = A(:,end)==0; %...
U(id,:) = NaN;
%         %quality control
%         fname = '../clean_results/bigvolumeIS_tr.raw';
%         ny = 1601;
%         nx = 1645;
%         nz = 400;
%         V = lect_raw(fname,nx,ny,nz,5);
%         figure;imshow3D(permute(V,[1 3 2]))
%         
%         tmp_U = reshape(U(:,3),siz_DVC)-mean(U(~isnan(U(:,3)),3));
%         figure;imshow3D(permute(tmp_U,[3 1 2]),[-1e-3 1e-3])
%         isl = 80;
%         tmp_U_posX = reshape(U_pos(:,1),siz_DVC);
%         tmp_U_posY = reshape(U_pos(:,2),siz_DVC);
%         tmp_U_posZ = reshape(U_pos(:,3),siz_DVC);
%         tmp_X = squeeze(tmp_U_posZ(:,isl,:));
%         tmp_Y = squeeze(tmp_U_posX(:,isl,:));
%         tmp_C = squeeze(tmp_U(:,isl,:));
%         figure;imagesc(tmp_X(:),tmp_Y(:),tmp_C);
%         axis equal;axis tight
%         id = abs(Gav(:,2)-isl*incDVC)<100;
%         figy = Gav(id,1);
%         figx = Gav(id,3);
%         hold on;plot(figx,figy,'.')
%         colorbar
        

%U = inpaint_nans(U);
%fill NaN - see Ahmet code
%Not a good idea, create some unusual high peak

U_xpos = reshape(U_pos(:,1),siz_DVC);
U_ypos = reshape(U_pos(:,2),siz_DVC);
U_zpos = reshape(U_pos(:,3),siz_DVC);
    %figure;imshow3D(U_zpos)

n = length(Nav(:,1));
m = length(U(:,1));

U_interp = zeros(n,2*lambda+1);
peak_shape = zeros(n,2*lambda);
peak_value = zeros(n,1);
inor = -lambda:1:lambda;
inor2 = -lambda:1:lambda-1;
figure
tic
%for j = 1:n
for j = 10000:10050
    U_nor = U(:,1).*Nav(j,1) + U(:,2).*Nav(j,2) + U(:,3).*Nav(j,3);
    U_nor = reshape(U_nor,siz_DVC);
    X_interp = inor.*Nav(j,1) + Gav(j,1);
    Y_interp = inor.*Nav(j,2) + Gav(j,2);
    Z_interp = inor.*Nav(j,3) + Gav(j,3);
    U_interp(j,:) = interpn(U_xpos, U_ypos, U_zpos, U_nor, X_interp, Y_interp, Z_interp, 'linear');
    hold on;plot(inor,U_interp(j,:))
    %derive the profile
    %peak_shape(j,:) = diff(U_interp(j,:),1,2);
    %hold on;plot(inor2,peak_shape(j,:))
    %stock the peak value
    %peak_value(j) = max(peak_shape(j,:));
end
toc