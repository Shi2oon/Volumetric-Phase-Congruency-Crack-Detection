%Convert the raw data into Davis-input format through post-processing
%15/05/2019

%% add paths
addpath fcts\

%% read image sequence
%dirname : data directory
%valeurs : selected slices range

dirname = '../../CN_01_450N_287kcycles/Reconstruction_cluster_Full';
valeurs = [851, 1250];
A = readImSeq(dirname,'tif',valeurs(1),valeurs(2));
siz = size(A);

%% convert float format to uint8 format

tmp = A(:,:,ceil(siz(3)*0.45):ceil(siz(3)*0.55));
    figure;hist(tmp(:),1000)

%vmin/max : see histogram

vmin = -1e-3;
vmax = 1.5e-3;
A8bit = ImFmtConvert(A(:,:,:),'uint8',[vmin,vmax]);
siz = size(A8bit);
    figure;imshow3D(A8bit)


%% manual rotation of an image
% ---- rotate in XY-plane (along Z axis)
        SL = A8bit(:,:,ceil(siz(3)/2));
            figure;imshow(SL)
        % pick two points
        [x,y] = ginput2(2);
            hold on;plot(x,y,'-og','LineWidth',2);
        % compute the rotation angle
        ang = atan( (y(2)-y(1))/(x(2)-x(1)) ) * 180/pi; %[degree]
        %save/load : keep the same angles for different images in the same 
        %            dataset
        save('ang-450N.mat','ang')
        load('ang-450N.mat')
% rotation procedure
SL = A8bit(:,:,1);
SL = imrotate(SL,ang,'bilinear','loose');
sizROT = size(SL);    sizROT = [sizROT(1), sizROT(2), siz(3)];
fprintf('rotating the image ...\n');  tic
Vrot = zeros(sizROT,class(A8bit));
for i=1:siz(3)
    SL = A8bit(:,:,i);
    SL = imrotate(SL,ang,'bilinear','loose');
    Vrot(:,:,i) = SL;
end
T=toc;fprintf(['image rotation finished, T=',num2str(T),' seconds\n']);
    figure;imshow3D(Vrot)


%% manual cropping
%Vrot : value must be found on the previous figure

A8bit_cropped = Vrot(651:2250,851:2450,:);
    figure;imshow3D(A8bit_cropped)
    
%% saving
%fnamesS : includes the repertory and the name of the file for saving

fnameS = '../Raw_Data/cast_iron/450N/X1600-Y1600-Z400-287kcycle.raw';
saveImRAW(fnameS,A8bit_cropped);
