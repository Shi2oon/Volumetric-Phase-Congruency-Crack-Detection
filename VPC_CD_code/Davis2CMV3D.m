%Convert the DVC result from Davis to the format compatible for CMV3D
%07/05/2019
%-Yang Chen


%% read DVC data

vxsiz = 1e-3; %voxel size in mm

%fname : includes the repertory and the name of the file for loading
fname = ['Results_DaVis\Pablo_CastIron\450N\50k-190k-250k-270k-287k_32laststep_75_overlap\287k.dat'];

%Ly : offset equal to the y-axis length
Ly = 1600;

A = dlmread(fname,' ',3,0);
nbr_pts = size(A,1);


%% kill bad correlation 
% id = A(:,7)==1;
% A = A(id(:),:);
% nbr_pts = size(A,1);
% nx = unique(A(:,1));
% nx = size(nx,1);
% ny = unique(A(:,2));
% ny = size(ny,1);
% nz = unique(A(:,3));
% nz = size(nz,1);


%% get the number of points

fid=fopen(fname,'r');
ligne = fgetl(fid);
ligne = fgetl(fid);
ligne = fgetl(fid);
fclose(fid);
nx = str2num(ligne(strfind(ligne,'I=')+2:strfind(ligne,', J=')));
ny = str2num(ligne(strfind(ligne,'J=')+2:strfind(ligne,', K=')));
nz = str2num(ligne(strfind(ligne,'K=')+2:strfind(ligne,', F=')));

%% reference grid

refX = A(:,1)./vxsiz;
refY = A(:,2)./vxsiz;
refZ = A(:,3)./vxsiz;

% Visual verification
%
%     inc = 10;
%     figure;plot3(refX(1:inc:end),refY(1:inc:end),refZ(1:inc:end),'.')
%     hold on;quiver3(refX(1:inc:end),refY(1:inc:end),refZ(1:inc:end),...
%                    A(1:inc:end,4),A(1:inc:end,5),A(1:inc:end,6),1000)
%     xlabel('x')
%     ylabel('y')

%% deformed grid

defX = refX + A(:,4)./vxsiz;
defY = refY + A(:,5)./vxsiz;
defZ = refZ + A(:,6)./vxsiz;
PH = A(:,end);

%% arrange the data into CMV3D-compatible format

A_ref = [[1:nbr_pts]',ones(nbr_pts,1),refX,-(Ly-refY),-refZ];
A_def = [[1:nbr_pts]',PH,defX,-(Ly-defY),-defZ];

% NOT NECESSARY
% re-number the data first following line (y), then colomn(x), then slice(z)
% --> in order to use "regular grid" in CMV3D
%     isort = reshape([1:nx*ny*nz],[nx,ny,nz]);
%     isort = permute(isort,[2 1 3]);
%     isort1 = reshape([1:nx*ny*nz],[ny,nx,nz]);
%     isort = isort(isort1);
% A_ref = A_ref(isort(:),:);
% A_ref(:,1) = [1:nbr_pts]';
% A_def = A_def(isort,:);
% A_def(:,1) = [1:nbr_pts]';
% A_ref = round(A_ref);

%% save to CMV3D-compatible format

fnameS = [fname(1:end-4),'_ref.pts'];
fileID = fopen(fnameS,'w');
fprintf(fileID,'%i\t %i\t %i\t %i\n',[nbr_pts,nz,ny,nx]);
fprintf(fileID,'%i\t %i\t %.12f\t %.12f\t %.12f\n',A_ref');
fclose(fileID);
fnameS = [fname(1:end-4),'_def.pts'];
fileID = fopen(fnameS,'w');
fprintf(fileID,'%i\t %i\t %i\t %i\n',[nbr_pts,nz,ny,nx]);
fprintf(fileID,'%i\t %i\t %.12f\t %.12f\t %.12f\n',A_def');
fclose(fileID);

