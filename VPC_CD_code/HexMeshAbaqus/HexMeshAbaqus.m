function [Nodes, Elements, StepBCs,offset, dcdataO, missingnodes,VertMskCloud, ...
     Faces_MaskPtCloud] = HexMeshAbaqus(dvcdata,BinaryMask_3D,INPfolder,unit)
close all;
 set(0,'defaultAxesFontSize',25);       set(0,'DefaultLineMarkerSize',14) 
% BinaryMask_3D is 1 and 0 mask, zero is for no masked area

[dvcdata, ~, offset] = unist_3D_Abaqus(dvcdata,unit);

% organise data and getting steps
X  = dvcdata(:,1);              x = unique(X);          Sx = abs(x(2)-x(1));
Y  = dvcdata(:,2);              y = unique(Y);          Sy = abs(y(2)-y(1));
Z  = dvcdata(:,3);              z = unique(Z);          Sz = abs(z(2)-z(1));
Ux = reshape(dvcdata(:,4),length(x), length(y),length(z));
Uy = reshape(dvcdata(:,5),length(x), length(y),length(z));
Uz = reshape(dvcdata(:,6),length(x), length(y),length(z));
Dis_3D(:,:,:,1) = Ux;
Dis_3D(:,:,:,2) = Uy;
Dis_3D(:,:,:,3) = Uz;
if isempty(BinaryMask_3D)
	[a, b, c, ~ ] = size(Dis_3D);
    BinaryMask_3D = ones(a, b, c);
end
% add Nans to mask
[a]=isnan(squeeze(Dis_3D(:,:,:,1)));
BinaryMask_3D(a)=0;
Plot3D(BinaryMask_3D,[],[],[],[],'3D Binary Mask');     
close;
% plotAllDis(dvcdata,unit,'Uz'); close
% Make point clouds of voxelized data
[BinaryMask_PtCloud, XYZmskDis] = MakePtClouds(BinaryMask_3D,Dis_3D);

% Create Hexahedral mesh in ABAQUS format from voxelized data
[VertMskCloud, Faces_MaskPtCloud] = HexMesh3D_Abaqus8Node(BinaryMask_PtCloud);

% Display Hexahedral mesh
close all; DisplayAbaqusMesh(VertMskCloud, Faces_MaskPtCloud); 
saveas(gcf,[INPfolder '\Hex_Mesh.tif'],'tiffn');  close   
% saveas(gcf,[INPfolder '\Hex_Mesh.fig']);  
dcdataO = [XYZmskDis(:,1).*Sx  XYZmskDis(:,2).*Sy  XYZmskDis(:,3).*Sz ...
           XYZmskDis(:,4)      XYZmskDis(:,5)      XYZmskDis(:,6)];
VertMskC = [VertMskCloud(:,1).*Sx  VertMskCloud(:,2).*Sy  VertMskCloud(:,3).*Sz];
       
% Write data to file .. need to modfiy the writting to resemble abaqus
[Nodes, Elements, StepBCs, missingnodes] = CreateINPfilefrom3DDataset ...
    (VertMskC, dcdataO, INPfolder);

save([INPfolder '\Meshing.mat']);
end