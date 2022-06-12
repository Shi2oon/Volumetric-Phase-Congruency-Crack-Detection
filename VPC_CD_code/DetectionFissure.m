%Crack segmentation from the subtracted image (output of CMV3D)
%Summer 2019

%% add path
addpath fcts

%% read image
%fname : includes the repertory and the name of the file for loading
fname = '../Results_CMV3D/cast_iron/450N/190k/190k_IS.raw';
%nx, ny, nz = width, height and number of slice
nx = 1600;
ny = 1600;
nz = 400;
A = lect_raw([fname(1:end-4),'.raw'],ny,nx,nz,5);
siz = size(A);
typeA = class(A);
    figure;imshow3D(A);
    
%% reduce the noise

hfilter = 4;    %parameters to be optimized
stdfilter = 1.5; %parameters to be optimized
Af = zeros(siz,typeA);
for i=1:siz(3)
    h = fspecial('gaussian',hfilter,stdfilter);
    Af(:,:,i) = imfilter(A(:,:,i),h);
end
    figure;imshow3D(Af)   
    
%% crack threshold

% for i=1:50
%     thres0 = 70+i; %parameter to be optimized
%     id = Af<thres0;
%     nbrpts = sum(id(:)==1);
%     Diff = (numel(Af) - nbrpts);
%     OptiTresh(i) = Diff;
%     threshold(i) = thres0;
% end
%     figure;plot(threshold,OptiTresh);
%     xlabel('Threshold Value')
%     ylabel('Number of points eliminated')
%     
% text = 'Choose finalthresh value according to this graph \n';     
% finalthresh = input(text); 

finalthresh = 100; %parameters to be optimized
id = A < finalthresh;
    close all
    figure;imshow3D(id)
    
% Visual verification
% isl = 223;
% [figy,figx] = find(id(:,:,isl));
% figure;imshow(A(:,:,isl))
% hold on;plot(figx,figy,'.')
   
%% manual mask
% only keep the values inside the rectangle created by the two selected
% points
[pylon1] = ginput(1);
[pylon2] = ginput(1);
fenetre_x = round([min(pylon1(2),pylon2(2)),max(pylon1(2),pylon2(2))]);
fenetre_y = round([min(pylon1(1),pylon2(1)),max(pylon1(1),pylon2(1))]);
maskedA = zeros(size(id));
maskedA(fenetre_x(1):fenetre_x(2),fenetre_y(1):fenetre_y(2),:) =...
     id(fenetre_x(1):fenetre_x(2),fenetre_y(1):fenetre_y(2),:);

     close all
     figure;imshow3D(maskedA)

%% particle analysis

conn = 26;
analysedA = bwareaopen(maskedA,40,conn); %parameters to be optimized
    figure;imshow3D(analysedA)
    
%  for i=1:30
%     objlen = i; %parameter to be optimized
%     analysedA = bwareaopen(maskedA,objlen,conn);
%     pointamo = sum(analysedA(:)==1);
%     gap = (numel(maskedA) - pointamo);
%     OptiPix(i) = gap;
% end
%     figure;plot(OptiPix);
%     xlabel('Number of pixels in an object')
%     ylabel('Number of points eliminated')
%     title('Choose P value according to this graph')   

%% kill all point far from the middle slice (crack location)

analysedA(:,:,1:50) = 0;
analysedA(:,:,350:end) = 0;
close all;
figure;imshow3D(analysedA)

%% save image
%fileName : includes the repertory and the name of the file for saving
fileName = '../Rdy_for_Avizo/cast_iron/450N_190k_rdy.raw';
saveImRAW(fileName,uint8(permute(analysedA,[2 1 3])));
 
%% quality control

close all
figure;
for i = 1:siz(3)
    imshow(A(:,:,i))
    [figy,figx] = find(analysedA(:,:,i)==1);
    hold on;plot(figx,figy,'.r')
    pause
end
