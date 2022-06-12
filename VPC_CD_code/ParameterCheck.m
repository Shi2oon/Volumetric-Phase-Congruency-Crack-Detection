clear all;

PC1 = 5;
PC2 = 1.5;
Med_Circle = 2;
Med_Thresh = 0.008;

inDir = 'X:\Shaocheng Ma\resliceddata\resliceonefor5Nvs37N_DVC_Vec_64x64x64_75_ov-30volumefractionand0.4\';
inFile1 = 'B00001.dat'; %%%% Data!

temp1=importdata([inDir inFile1]);
D = temp1.data(:,1:6);

[X1, Y1, Z1, Ux1, Uy1, Uz1] = V2D(D);


Ux1(Ux1==0) = nan;
Uy1(Uy1==0) = nan;
Uz1(Uz1==0) = nan;

DisplacementField1 = Uz1(:,:,:); % select larges magnitude signal

xaxis = 1:size(DisplacementField1,1);
yaxis = 1:size(DisplacementField1,2);
zaxis = 1:size(DisplacementField1,3);


%load('test.mat')
%%


for i = yaxis
    
    M1 = squeeze(DisplacementField1(:,i,:)); %%% Extract itterative slice, Cropped Volume
    M2(:,i,:) =inpaint_nans(M1,1); % fill in missing data with nans
    
    M_Del(:,i,:) = OutDel(squeeze(M2(:,i,:)),Med_Circle,Med_Thresh); %20, 0.001 % outlier deletion
    
    M(:,i,:) =inpaint_nans(squeeze(M_Del(:,i,:)),1); % fill in outlier deleted data
    
    VPC1(:,i,:) = phasecongmono(squeeze(M(:,i,:)), 6, PC1, 2, 0.55, 1.0,0.25,10,PC2,-1); % use phase congruency
    
    
end

i_vis = 55;
subplot(2,3,1)
imagesc(squeeze(DisplacementField1(:,i_vis,:))); axis image;
subplot(2,3,2)
imagesc(squeeze(M2(:,i_vis,:))); axis image;
subplot(2,3,3)
imagesc(squeeze(M_Del(:,i_vis,:))); axis image;
subplot(2,3,4)
imagesc(squeeze(M(:,i_vis,:))); axis image;
subplot(2,3,5)
imagesc(squeeze(VPC1(:,i_vis,:))); axis image;


VPC1(isnan(DisplacementField1))=0;
subplot(2,3,6)
imagesc(squeeze(VPC1(:,i_vis,:))); axis image;

%[z,x] = ginput(1);
%centre_point = [round(x) i_vis round(z)]