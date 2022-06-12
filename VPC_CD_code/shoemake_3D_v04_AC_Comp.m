function [Euler, DVC_data_new] =shoemake_3D_v04_AC_Comp(inDir,inFile)

%% shoemake_3D_v04_07 (alpha)
%{
%     Removes 3-D rigid body translation and rotation from a volumetric
%     dataset of displacement data (input as DAT-type format, as from
%     DaVis). In addition to data outputs for later analysis, several
%     visualisation options are available, and may be selected in Section 2
%     ("User Inputs"), including options to realign all the data to an
%     arbitrary reference frame, export the data as a DAT-type file,plot
%     orthogonal slices, bin data and plot as separate orthogonal slices
%     and plot variation in average and standard deviation in z-direction.
%     The majority of inputs for code may be entered in "User inputs"
%     section, but some modifications must be made within the code.
%
% INPUT DATA FORMAT (DAT-type)
%   Data should be in a 6 column matrix of
%        [Coordinates(X,Y,Z) Displacements(dX,dY,dZ)]
%   Coordinates should be ordered to agree with:
%       [YY XX ZZ] = meshgrid(Ymax:-1:Ymin, Xmin:Xmax, Zmin:Zmax)
%       Coordinates = [XX(:) YY(:) ZZ(:)]  
%
% OUTPUTS
%   Structures: initial     {Data, Size, DX}
%               C           {Data, Size, DX}
%               DVC         {Data, Size, DX}
%               strainAll   {Data, Size, DX}
%               bin         {Data, Size, DX}
%               binStrain   {Data, Size, DX, Comp}
%
% Code originally developed by M. Mostafavi, with additions by M. Jordan.
% Based on Ken Shoemake's Euler angle extraction.
%
% Last edit: April 2014 (M Jordan)
% (C) 2014, University of Oxford
%}

%% 1. Routine Stuff
% format compact
% clc 
% clear
% set(0,'DefaultFigureWindowStyle','docked')
% close all; 
% pause off
version = 'shoemake_3D_v04.07 (alpha)';
disp(['Running ' version '.']); 
disp ' ';

%% 2. (a) User inputs: details of dataset and outputs
%Path for input file (.dat)
% inDir = ['U:\DVC_Data\1.DVC_UD_Forward_NoMica\Raw_x0_5\DVC_Vec_128x128x128_50%ov\'];

% inFile = 'D1.dat';
% cycleName = 'test';

% Overwrite variable "data" already in memory and repeat calculations? 
Overwrite = true;          %(true/false)

% Dataset information not recoverable from input file
res = 1.8e-3;               % resolution in mm per pixel
ov = 50;                    % The percentage overlap setting in final pass

% Crop notch? (For EE8519 only) 
CropNotch = false;           %(true/false)
notchDiameter = 4;          % 4mm at -8.1d, 2mm at -6.6d

% Crop borders (removes spurious DVC values at boundaries)
CropBorders = false;         %(true/false)
borderSize = 2;             %Size (in voxel) of border effects

% Correct for rigid body motions
CorrectRBM = true;          %(true/false)

% Reference displacements to geometry? 
GeometryReference = true;   %(true/false)
forceTheta = [30 10 20];   %**NI** 3 angles along principal axes in degrees e.g. [45 -45 60] (default:0 or false)

% Bin data? 
BinData = false;            %(true/false)
numberOfBins = 5;
PlotBinnedSlice = true;        %(true/false) Requires PlotOrthoSlices

% Extract orthogonal strains?
ExtractStrain = false;      %(true/false)

% Output options (true/false)
WriteDataToFile = false;       %(true/false)
outFolder = 'Test_data';
writeRotated = true;
printHeader = false;
requestWriteAll = true;       % (true/false) as a .dat type file, as from DaVis
requestedComponents = [];   % 3D coords + {Ux,Uy,Uz = 1,2,3 respectively}

PlotOrthoSlices = true;       %(true/false) If TRUE go to 2.1

PlotAverageAndStDev = false;   %Calc & plot the x-y slice average variation     
%>>Exclude: [#slicesFromTop #slicesFrombottom]: [0 0] is entire data set
exclusionRange = [0 0];
xScale = 'mm';                     %x-axis in 'pixel' or 'mm'

%% 2.1 Orthoslice options
% List requested plots in 4 element row vectors: 
%{
% PLOTS = {[orientation,sliceNo,map,Vectors] 'Title'}
% ORIENTATION (for map): {Slice i=sliceNo, where i=1-3 (X-Z)};
% SLICENO: valid range: 1 <= sliceNo <= lx,ly,lz (depends on orient);
Fractional values (0<=sliceNo<1) will also work for the whole volume.
% MAP: {1-3 = Ux-Uz, 4-6 = Exx-Ezz} N.B. 4-6 will require manual alteration
% of code (search 'mapData' and change variable name as appropriate);
% Vectors: true/false; 
% Placeholder 'Title' generates either Default figure title or 'string' as required;
%
% e.g. Plots = {...
%     [2 100 3 true] 'Title';...
%     [1 50 2 false] 'A';...
%     };
% generates 2 plots: 
% 1st: XZ slice (Y=100), Uz map, vectors & default title, 
% 2nd: YZ (X=50),Uy,no vectors and titled 'A' 
%}

% %Either specify particular plots ('Plots' can be any number of plots >0)
Plots = {...
    [3 0.5 2 true] 'Title';...
    [2 0.5 2 true] 'Title';...
    [1 0.5 2 true] 'Title';...
% %     [3 19 1 true] 'Title';...
%     [3 0.5 2 true] 'Title';...
%     [3 0.5 1 true] 'Title';... 
%     [2 0.5 2 true] 'Title';...
%     [1 0.5 2 true] 'Title';...
    };

% %Or for a series of plots, uncomment the following:
% first = 2; int = 1; last = 4; 
% plot_nos = first:int:last;
% Plots = repmat({[1 0 1 true] 'Title'}, size(plot_nos,2),1);
% for i= 1:size(plot_nos,2)
%     Plots{i,1}(2) = plot_nos(i);
% end

% Plot separate figure with original input data for comparison?
PlotOriginal = true;  %(true/false)

% Specify an offset origin for the plots 
%Coordinate origin options: 'Natural', 'Centre', 'Rotation' or arbitrary [x y z] in mm;
coordinateOrigin = 'Centre';
%Displacement origin options: 'Absolute', 'Relative' or arbitrary [dx dy dz] in mm;
displacementsOrigin = 'Relative';

%For map
% Smoothing function
% F = 1; 
F = fspecial('gaussian', 3, 0.5);
% Number of contours
CN = 24;
%For vector field
% Plot every nth vector on a square array or specify as [nx ny nz]
n =5;
% Arrow head size
A_S = 1;
% Vector colour
A_C = [0 0 0];
% Arrow thickness
A_T = 2;
% Arrow size
V_S = 2;
% Units for axes
units = 'mm';   %'mm' or 'voxel'

outOpts = {F,CN,n,A_S,A_C,A_T,V_S,units};

%% End of Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3. Import data, clean and remove rigid body motion
% Data import and calculations are slow, and so will only occur if new data
% is unique or overwrite requested. END is currently just before geometry
% reference section.

if ~exist('data','var') || Overwrite || ~strcmp(datPath,[inDir '\' inFile]);
    %% 3.1 Import data (if required) and data size
    datPath = [inDir '\' inFile]
    
%     inDir
    cd(inDir);
%     disp(['Importing data for ' cycleName])
    importDat(datPath);

%     % Extract size of data from heading line of .dat file
%     %> if error thrown comment out and uncomment alternative method
%     lx = str2double(colheaders{1,4}(3:end-1));    %NB curly brackets      
%     ly = str2double(colheaders{1,5}(3:end-1));    %NB curly brackets
%     lz = str2double(colheaders{1,6}(3:end));      %NB, runs to 'end' not 'end-1'
%%
D = data;    
j=-1;
i=1;
while j<0
   if D(i,2)==D(i+1,2)
       i=i+1;
   else
       lx=i;
       j=1;
   end
end
lx;

tempZ=reshape(D(:,3),lx,size(D,1)/lx);
%%

j=-1;
i=1;
while j<0
   if tempZ(1,i)==tempZ(1,i+1)
       i=i+1;
   else
       ly=i;
       j=1;
   end
end
ly;
lz=size(D,1)/ly/lx;
%%

    % Alternative method:extract size of the data when no header line
%     [lx ly lz dims] = bruteSize(tdata);
            
    % Distance between interrogation points in pixels
    dx_mm = data(2,1)-data(1,1);
    dx = dx_mm/res;
    % dx = (100/(100-ov))*(data(2,1)-data(1,1))/res*10^3;
    
    % Collect current data, dataset size and voxel size together
    %> To access information stored in a struct, use e.g. C.Data or C.Size
    C = struct('Data',data,'Size',[lx ly lz],'DX', dx_mm);   %C for current
    disp(['Size of data /voxel: (lx ly lz) = (' num2str(C.Size) ').']);    
    
    % Definitions
    component = {'U_x','U_y','U_z','|U|'};
    outDir = [pwd outFolder];
    eulerAngles = [nan;nan;nan];
    notchOrigin = [0,0];
    theta_deg = [0 0 0];
    
    %% 3.2 Crop notch (only for EE8519)


%     
    %% 3.3 Cleaning: Replace rows with 0 displacements to NaN
    
    tdata_x = C.Data(:,4);
    tdata_y = C.Data(:,5);
    tdata_z = C.Data(:,6);
    notn = (tdata_x==0 & tdata_y==0 & tdata_z==0);
    tdata_x(notn) = NaN;
    tdata_y(notn) = NaN;
    tdata_z(notn) = NaN;

    C.Data(:,4:6) = [tdata_x tdata_y tdata_z];
    
    %% 3.3.1 Removing border effects
    if CropBorders
        smallData = C.Data; 
        smallData = smallData(smallData(:,1)> borderSize*C.DX & smallData(:,1)<(C.Size(1)-borderSize)*C.DX,:);
        smallData = smallData(smallData(:,2)> borderSize*C.DX & smallData(:,2)<(C.Size(2)-borderSize)*C.DX,:);
        smallData = smallData(smallData(:,3)> borderSize*C.DX & smallData(:,3)<(C.Size(3)-borderSize)*C.DX,:);

        C.Data = smallData;
        C.Size = C.Size-2*borderSize;
    end
    
    initial = C;
    %% 3.3.2 Clean dataset, and correct for rigid body motions
    if  CorrectRBM
        [C.Data,eulerAngles,rotCentre] = RBMCorr('Minimum',C.Data,C.Size);
    else 
        disp 'No corrections for RBM made.'
    end
        
    %{    
    %% 3.3.2 Create a clean data (cdata) by deleting the rows with NaN
    cdata=C.Data(~any(isnan(C.Data),2),:);

    if(sum(sum(sum(isnan(cdata))))) ~=0
        disp('WARNING: data contains vectors with partial NaN entries');
    end

    %% 3.4 Remove rigid body motions 
    %% 3.4.1 RBD from displacements
    C0=zeros(1,6);
    C0(1,4:6)=nanmean(cdata(:,4:6));      %row of 3 RBD coordinates
    U0=repmat(C0,size(cdata,1),1);
    cdata=cdata-U0;

    %% 3.4.2 Finding the centre of rotation and setting as coordinate origin
    % Defining as the row with the smallest displacement after subtracting the
    % rigid body motion and setting it as the origin

    [X0_disp,X0L]=min(sum(abs(cdata(:,4:6)')));     %minium vector of cdata
    X0=zeros(1,6);
    X0(1,1:3)=cdata(X0L,1:3);                       %X0 coords = index coords of rotn centre
%     X0(1,1:3)=[5 2 2];
    X0=repmat(X0,size(cdata,1),1);
    cdata=cdata-X0;

    %xi0 is the location of points before loading and xip is the location after
    %displacement

    xi0=cdata(:,1:3);                   % Reference coords
    xip=cdata(:,1:3)+cdata(1:(end),4:6);      % Comparison dataset coords

    %% 3.4.3 Calculating rotation matrix, and checking the rank 
    % Rank should equal 3 to avoid singularity
    rank_xi0 = rank(xi0');
    rank_xip = rank(xip');

    if rank_xi0 ~= 3
        disp(['xi0 rank incorrect. rank_xi0 = ' num2str(rank_xi0)])
        ERROR
    end
    if rank_xip ~= 3
        disp(['xip rank incorrect. rank_xip = ' num2str(rank_xip)])
        ERROR
    end

    m = xi0\xip;  % backslash is MATLAB matrix division - solves for rotation efficiently 
    
    %% 3.4.4 Extracting Euler angles from the rotation matrix
    t=zeros(3,1);

    t(1)=atan2(m(2,3),m(3,3));
    cos_t2=sqrt(m(1,1)^2+m(1,2)^2);
    t(2)=atan2(-m(1,3),cos_t2);
    t(3)=atan2(m(1,2),m(1,1));

    disp(['Euler angles /degrees: (alpha beta gamma) = (' num2str(radtodeg(t'),'% 1.3f ') ').'])

    %% 3.4.5 Subtract displacements due to the pure (theoretical) rotation
    disp '##RB rotation correction disabled##'

    c=cos(t);
    s=sin(t);

    % theoretical rotation matrix based on the calculated angles

    m_theo=[c(2)*c(3), c(2)*s(3), -1*s(2); ...
         s(1)*s(2)*c(3)-c(1)*s(3), s(1)*s(2)*s(3)+c(1)*c(3), s(1)*c(2); ...
         c(1)*s(2)*c(3)+s(1)*s(3), c(1)*s(2)*s(3)-s(1)*c(3),c(1)*c(2)];

     xip_theo=(m_theo*xi0')';
     disp_theo=-xip_theo+xi0;       %= orig coords - rotn coords

    %% 3.4.6 Recreate original and cleaned displacement fields
    % Ux,Uy and Uz are displacements. Parameters without a prefix are
    % original and rb means rotated back. big_new_data_natural is the final
    % matrix similar to data corrected for rotation referenced to input x-y
    % frame.

    % rUi = disp_theo(:,1:3)  {theoretical displacments} 
    % Ui = cdata(:,4:6)       {measured displacements}
    % rbUi = Ui - rUi         {deformation displacements}
    
    rbUi = cdata(:,4:6) - disp_theo(:,1:3); 
    new_data=[cdata(:,1:3) rbUi];      %[CleanedCoords DeformDisps] N.B. cdata has been recentred c.f. tdata

    big_new_data_natural=C.Data;

    j=1;
    for i=1:size(big_new_data_natural)
        if ~isnan(big_new_data_natural(i,6))
            big_new_data_natural(i,:)=new_data(j,:);
            j=j+1;    
        end
    end
    %}
    %% 3.4.7 Reform to DaVis-like file (inc 0 displacement coordinates)
    if CropBorders
        C.Data(:,1:3) = smallData(:,1:3);
    else
        C.Data(:,1:3) = data(:,1:3);
    end

    DVC_new_data=C.Data;
    DVC_new_data(isnan(DVC_new_data))=0;
    
    DVC = struct('Data',DVC_new_data,'Size',C.Size,'DX', dx_mm);
    
    correctedData = C.Data;
    disp(['Size of corrected data (pre-rotation) /voxel (lx ly lz) = (' num2str(C.Size) ').']);
else
    disp('*** Using previous results until geometry reference section (S4.1) ***')
    C.Size = initial.Size;
    C.DX = initial.DX;
    disp('Dataset details:')
    disp(['   Size of corrected data /voxel: (lx ly lz) = (' num2str(C.Size) ').']);    
    disp(['   Euler angles /degrees: (alpha beta gamma) = (' num2str(radtodeg(eulerAngles),'% 1.3f ') ').'])
    disp('****')
end
disp(' ')

%% 4 Rotate coordinates and displacements to align with geometry
C.Data = correctedData;
% theta_deg=[0 30 0]; disp(['theta_deg overwritten (= ' num2str(theta_deg) ')'])
theta_deg = eulerAngles
if GeometryReference && (mod(sum(abs(theta_deg)),360) ~= 0);
%     theta_deg
    [C.Data, C.Size, C.DX] = imRotateShell(C.Data, C.Size, C.DX, theta_deg);
end

DVC_data_new = C;
Euler = eulerAngles; 
    %{
if GeometryReference && (mod(sum(abs(theta_deg)),360) ~= 0)    
    
    [C.Data(:,4:6)] = dispNotchRef(C.Data(:,4:6), theta_deg);
% OC6 % OC8
    [C.Data, C.Size, C.DX] = imRotateShell(C.Data,C.Size,C.DX,theta_deg);
    end
%}

% %% 5. Quantification
% %% 5.1 Strain calculation
% if ExtractStrain
%     %An inelegant method of collecting Exx,Eyy and Ezz and coordinates in a
%     %columnwise data format. calcOrthoStrain is embedded in plotOrthoSlice
%     %(S.7) and so it is not necessary to compute the strains here if only
%     %wish to plot them. If binning data, strains calculated here may be binned afterwards(S5.2),
%     %while strains calculated in plotOrthoSlice will be from binned data.
%     
%     [dataOut,coordsOut,sizeOut,compOut] = calcOrthoStrain(4,C.Data,C.Size,'simple',C.DX);
%     strainDataOut = zeros(prod(sizeOut),6);
%     strainDataOut(:,1:4) = [coordsOut dataOut];
%     
%     for i=5:6
%         [strainDataOut(:,i),coordsOut,sizeOut,compOut] = calcOrthoStrain(i,C.Data,C.Size,'simple',C.DX);
%     end
%     strainAll = struct('Data',strainDataOut,'Size',sizeOut,'DX',C.DX,'Component',['All']);    
% end
% % OC5
% %% 5.2 Data binning along notch direction
% 
% clear bin binStrain
% if BinData
% % OC7    
%     binDir = 1;
%     [dataOut,sizeOut,dxOut] = xbinData(binDir,C.Data,C.Size,C.DX,numberOfBins);
%     bin = struct('Data',dataOut,'Size',sizeOut,'DX', dxOut);
%     disp(['Binned data available >> using ' num2str(numberOfBins) ' bins.'])
%     
%     if ExtractStrain
%         binDir = 1;
%         [dataOut,sizeOut,dxOut] = xbinData(binDir,strainAll.Data,strainAll.Size,strainAll.DX,numberOfBins);
%         binStrain = struct('Data',dataOut,'Size',sizeOut,'DX', dxOut);
%         disp(['Binned strain data available >> using ' num2str(numberOfBins) ' bins.'])              
%     end
% else disp('No binned data available')
% end
% 
% %% 6. Write rotated data to file
% 
% if WriteDataToFile
%     if writeRotated
%         writeThis = C.Data;    
%         writeSize = C.Size;
%     else
%         writeThis = DVC.Data;
%         writeSize = DVC.Size;
%     end
%     
%     writeThis(isnan(writeThis))=0;    
%     writeDataToFile(requestedComponents,requestWriteAll,printHeader,outFolder,writeThis,writeSize,cycleName,eulerAngles)
%     clear writeThis writeSize
%     % OC1
% else disp 'Write data to file not requested.'
% end
% 
% disp(' ')
% 
% %% 7. Plotting of 2D ortho slices
% % OC2
% if PlotOrthoSlices
% % OC4
%     if strcmp(coordinateOrigin,'Rotation')
%         if exist('rotCentre','var')
%             coordOrigin = rotCentre;
%         else
%             disp 'Rotation centre underfined. Orthoslices centred on (0,0,0)'
%             coordOrigin = 'Natural';
%         end
%     else 
%         coordOrigin = coordinateOrigin;
%     end
%     
%     borderWidth = borderSize*CropBorders;
%     
%     for i=1:size(Plots,1)
%         %Slice orientation, number, requested map and vector option
%         mapID = Plots{i,1};  %[2 Y Map Vectors]; 
%         if strcmp(Plots{i,2},'Title')
%             plotTitle = {true;'(rot. corr.)'; cycleName};
%         else
%             plotTitle = {false;Plots{i,2}; ''};
%         end
%         disp(['Plot ' num2str(i) ': {[' num2str(mapID) '] ' Plots{i,2} '}...']);
%         
%         %Requested slice number is adjusted to account for rotation of
%         %sample
%         if mapID(2)<1
%             mapID_temp=ceil(mapID(2)*C.Size(mapID(1)));
%         else
%             mapID_temp=ceil(C.Size(mapID(1))*mapID(2)/initial.Size(mapID(1)));
%         end
% 
%         mapID_corr = [mapID(1) mapID_temp-borderWidth mapID(3:end)];
%         disp(['** Plotting corrected slice: ' num2str(mapID_temp) ' **'])
%         
%         figure
%         plotOrthoSlice(mapID_corr,C.Data,C.Size,...
%             coordOrigin,displacementsOrigin,outOpts,plotTitle,C.DX);
% %         caxis([-0.0084 0.0091]); disp '##caxis altered'
%         if PlotOriginal
%             % Plot original data for comparison
%             plotTitle = {true;'(original)'; cycleName};
% 
%             %Requested slice number is adjusted to account for rotation of
%             %sample
%             if mapID(2)<1
%                 mapID_temp=ceil(mapID(2)*initial.Size(mapID(1)));
%             else
%                 mapID_temp=mapID(2);
%             end
%             
%             mapID_init = [mapID(1) mapID_temp-borderWidth mapID(3:end)];
%             disp(['** Plotting original slice: ' num2str(mapID_temp) ' **'])            
%             
%             figure
%             plotOrthoSlice(mapID_init,initial.Data,initial.Size,...
%                 coordOrigin,displacementsOrigin,outOpts,plotTitle,initial.DX);
% %             caxis([0.0084 -0.0091]); disp '##caxis altered'
%         end
%         
%         if PlotBinnedSlice && exist('bin','var')
%             % Plot binned data
%             plotTitle = {true;'(binned)';cycleName};
%             outOptsBin = {F,CN,[1 n n],A_S,A_C,A_T,V_S,units};
%             %Requested slice number is adjusted to account for rotation of
%             %sample
%             if mapID(2)<1
%                 mapID_temp=ceil(mapID(2)*bin.Size(mapID(1)));
%             else
%                 mapID_temp=ceil(bin.Size(mapID(1))*mapID(2)/initial.Size(mapID(1)));
%             end
%             
%             mapID_bin = [mapID(1) mapID_temp-borderWidth mapID(3:end)];
%             disp(['** Plotting binned slice: ' num2str(mapID_temp) ' **'])
%             figure
%             plotOrthoSlice(mapID_bin,bin.Data,bin.Size,...
%                 coordOrigin,displacementsOrigin,outOptsBin,plotTitle,bin.DX);
% %             caxis([-0.0084 0.0091]); disp '##caxis altered'
%             % OC9
%         end
%     end
% else disp 'Ortho slices not requested.'
% end
% 
% disp(' ')
% %% 8. Plotting Average and StDev (slice-wise) in Z-direction
% 
% if PlotAverageAndStDev 
%     z_range = [exclusionRange(1)+1 C.Size(3)-exclusionRange(2)];
%     z_length = z_range(2)-z_range(1)+1;
% 
%     % RMS and Ave have column headings: 
%     %
%     % ||Slice_number|Z_coord|Blank!|<Ux>|<Uy>|<Uz>|<U>|{SD(U)/<U>|}|
% 
%     [RMS, ave] = stdevAndAverage(C.Data(:,4:6),C.Size,z_range,C.DX);   
%     %OC3
%     figure;  
%     hold on
%     switch xScale
%         case 'pixel'
%             xlabel('Distance along Z /pixel');
%             xlim(z_range);
%             xCol = 1;
%         case 'mm'
%             xlabel('Distance along Z /mm');
%             xlim(z_range*dx_mm);
%             xCol = 2;
%         otherwise
%             disp('ERROR: Invalid X-axis unit; valid units are {pixel/mm}')
%             ERROR
%     end
%     %     plot(Ave(:,xScale),Ave(:,4),'rx-',Ave(:,xScale),Ave(:,5),'gx-',Ave(:,xScale),Ave(:,6),'bx-',Ave(:,1),Ave(:,7),'k.-',[0,max(RMS(:,1))],[0,0],'k-');
%     errorbar(ave(:,xCol),ave(:,4),RMS(:,4)./2,'rx-');
%     errorbar(ave(:,xCol),ave(:,5),RMS(:,5)./2,'gx-');
%     errorbar(ave(:,xCol),ave(:,6),RMS(:,6)./2,'bx-');
%     errorbar(ave(:,xCol),ave(:,7),RMS(:,7)./2,'k.-');
%     plot([0,max(RMS(:,1))],[0,0],'k-');          %plot x-axis line at y=0
% 
%     title(['Average(U_i) over X-Y slices for ' cycleName]);
%     ylabel('X-Y Slice Average /mm');
%     legend('U_x','U_y','U_z','|U|','Location','EastOutside');
%     hold off
%     
% else disp 'Average and StDev plots not requested.'
% end
% %% 9. Save images (not implemented)
% 
% %{
% % save as tiff
% % whichpara=orig_Ez;
% % imheader='y:\tiffs\al-sic_orig_Ez';
% % num_zeros=3;
% % % create a hdf5 file for storing the sinograms
% % 
% % 
% % % Create the datasets
% % max_im=nanmax(whichpara(:));
% % 
% % for i=1:lz
% %     sprint_im=['%0' int2str(num_zeros) 'd'];
% %     imagename =[imheader sprintf('%04d',i) '.tif'];
% %     im=squeeze(whichpara(:,:,i))/max_im;
% %     imwrite(im,imagename,'tif','Compression','none');
% % end
% %}
% 
% %% 10. Cleaning up
% %{
% use >>CLEARVARS - EXCEPT VAR1 VAR2
% 
% clear A_C A_S A_T C C2 F U0 Ux Uy Uz V_S X X0L X0_disp Y c cUx0 cUy0 ...
%     cUz0 cX0 cY0 cZ0 cdata cleaned_Ez ... %cleaned_Ux cleaned_Uy cleaned_Uz ...
%     cleaned_x cleaned_y cleaned_z disp_theo h1 i j lx ly lz m m_theo n ...
%     oUx0 oUy0 oUz0 oX0 oY0 oZ0 orig_Ez  orig_x ...
%     orig_y orig_z rUx rUy rUz rbUx rbUy rbUz rbx rby rbz s sEz sUx ...
%     sUy sX sY sZ tdata tempEz tempUz1 tempUz2 tempZ temp_Ez x xi0 ...
%     xip xip_theo xx y yy z zz new_data big_new_data %orig_Uz
% %}
% 
% disp(' ')
% disp(['*** ' version ' complete ***'])
% % end
% shg
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Obselete Code (OC)
% % Indications of OC# relate to code removed from the main script and
% % generally replaced by a separate function. Kept for records until beta
% % testing past.
% 
% %{
% OC1 - replaced by separate function: writeDataToFile
%     if ~exist(outDir,'dir')
%         mkdir(outDir);
%         disp(['Saving rotated data in ' outDir])
%     else
%         disp('Output folder already exists; data will be overwritten.')
%         disp('=>Hit any key to continue')
%         pause
%         disp(['Saving rotated data in ' outDir])
%     end
%     
%     c=clock;
%     
%     % Write DVC_new_data
%     if requestWriteAll
%         disp('Writing DVC-all');
%         filename = [cycleName '-All.txt'];
%         outPath = [outDir filename];
% 
%         string1 = [cycleName ' - All rotated data calculated at ' num2str(c(4)) ':' sprintf('%02d',c(5)) ' on ' num2str(c(3)) '/' num2str(c(2)) '/' num2str(c(1))];
%         string2 = ['(I,J,K)=(' sprintf('%03d',dataSize(1)) ',' sprintf('%03d',dataSize(2)) ',' sprintf('%03d',dataSize(3)) ')']; 
%         string3 = ['|        x/mm|         y/mm|         z/mm|       ' component{1} '/mm|       ' component{2} '/mm|       ' component{3} '/mm|'];
% 
%         fileID = fopen(outPath,'w');
%         fprintf(fileID,'%s\r\n',string1);
%         fprintf(fileID,'%s\r\n',string2);
%         fprintf(fileID,'%s\r\n',string3);
%         fprintf(fileID,'%11.6E %11.6E %11.6E %11.6E %11.6E %11.6E \r\n',DVC_new_data');
%         fclose(fileID);
%     end
%     
%     if sum(requestedComponents)>0
%         disp('Writing displacement components to file');
%         % For compatibility with "Displacement_difference.exe" (Delphi file)
%         % require 4 column format:
%         %
%         % string1//
%         % string2//
%         % string3//
%         % x|y|z|Ui
% 
%         writeData = zeros(size(DVC_new_data,1),4);
%         writeData(:,1:3) = DVC_new_data(:,1:3);
% 
%         for i = requestedComponents   % change to 1:3 for all of x,y,z
%             filename = [cycleName '-' component{i} '.txt'];
%             outPath = [outDir filename];
% 
%             string1 = [cycleName '-' component{i} ' rotated data calculated at ' num2str(c(4)) ':' sprintf('%02d',c(5)) ' on ' num2str(c(3)) '/' num2str(c(2)) '/' num2str(c(1))];
%             string2 = ['(I,J,K)=(' sprintf('%03d',dataSize(1)) ',' sprintf('%03d',ly) ',' sprintf('%03d',lz) ')']; 
%             string3 = ['|        x/mm|         y/mm|         z/mm|       ' component{i} '/mm|'];
%             writeData(:,1:3) = DVC_new_data(:,1:3);
%             writeData(:,4) = DVC_new_data(:,i+3);
% 
%             fileID = fopen(outPath,'w');
%             fprintf(fileID,'%s\r\n',string1);       %probably change to one line
%             fprintf(fileID,'%s\r\n',string2);
%             fprintf(fileID,'%s\r\n',string3);
%             fprintf(fileID,'%11.6E %11.6E %11.6E %11.6E \r\n',writeData');
%             fclose(fileID);
%         end
%     end
% %%%%%%%%%%
% OC2 - This OC is not required due to 'Overwrite' if-loop which removes all need for calculation.
% 
% % % To select plots - if activated, overwrites Plots above.
% % % Either specify plots:
% % disp('>> Plots in section 2.b overwritten'); Plots={...
% %     [3 20 2 true] 'Title';...
% %     [ 100 2 true] 'Title';...
% %     [3 18 2 true] 'Title';...
% %     };
% 
% % % Or generate series of plots:
% % first = 1; int = 50; last = 213;
% % plot_nos = first:int:last;
% % disp('>>Plots in section 2.b overwritten')
% % Plots = repmat({[2 0 2 true] 'Title'}, size(plot_nos,2),1);
% % for i= 1:size(plot_nos,2)
% %     Plots{i,1}(2) = plot_nos(i);
% % end
% %%%%%%%%%%
% OC3 - previously required figures
%    
%     figure
%     plot(RMS(:,1),RMS(:,4)./Ave(:,4),'rx-',RMS(:,1),RMS(:,5)./Ave(:,5),'gx-',RMS(:,1),abs(RMS(:,6)./Ave(:,6)),'bx-',RMS(:,1),RMS(:,8),'k.-');
%     plot(RMS(:,1),RMS(:,4)./Ave(:,4),'rx-',RMS(:,1),RMS(:,5)./Ave(:,5),'gx-',RMS(:,1),abs(RMS(:,6)./Ave(:,6)),'bx-',[0,max(RMS(:,1))],[0,0],'k-');
%     title(['StDev/Ave(U_i) over X-Y slices for ' cycleName]);
%     xlabel('Distance along Z /mm');
%     ylabel('X-Y Slice (StDev/Average) /-');
%     legend('U_x','U_y','U_z','|U|','Location','EastOutside');
%     xlim(z_range*dx*res)
% 
% 
%     figure
%     Xfit = polyfit(Ave(1:55,4),RMS(1:55,4),1)
%     Yfit = polyfit(Ave(1:55,5),RMS(1:55,5),1)
%     Zfit = polyfit(abs(Ave(1:55,6)),RMS(1:55,6),1)
% 
%     plot(Ave(:,4),RMS(:,4),'rx',Ave(:,5),RMS(:,5),'gx',abs(Ave(:,6)),RMS(:,6),'bx');    %,[0,max(RMS(:,1))],[0,0],'k-');
%     hold on
%     plot([0,1e-3],[Xfit(2), Xfit(2)+Xfit(1)*1e-3],'r-',[0,1e-3],[Yfit(2), Yfit(2)+Yfit(1)*1e-3],'g-',[0,1e-3],[Zfit(2), Zfit(2)+Zfit(1)*1e-3],'b-')
%     title(['StDev vs Ave over X-Y slices']);
%     xlabel('|X-Y Slice Average| /mm');
%     ylabel('X-Y Slice StDev /mm');
%     legend('U_x','U_y','U_z','Location','EastOutside');
%     xlim([0, 1e-3]);
%     ylim([0, 1e-3]);
%     hold off
% 
%     figure
%     plot(RMS(:,1),RMS(:,4),'rx-',RMS(:,1),RMS(:,5),'gx-',RMS(:,1),RMS(:,6),'bx-',RMS(:,1),RMS(:,7),'k.-');
%     title(['StDev(U_i) over X-Y slices for ' cycleName]);
%     xlabel('Distance along Z /mm');
%     ylabel('X-Y Slice Standard Deviation /mm');
%     legend('U_x','U_y','U_z','|U|','Location','EastOutside');
%     xlim(z_range*dx*res)
% %%%%%%%%%%
% OC4 - for use with previous version of plotOrthoSlice
%     % Coordinates for cleaned map (rigid body displacement corrected)
% %     RBDdata = data(:,1:3);
% %     RBD = repmat(X0(1,1:3),length(data),1);
% %     RBDdata = RBDdata - RBD;
% %     mapData = cleaned_Ez;
% %     vectorData=big_new_data_natural(:,4:6);
% %     vectorData=big_new_data(:,4:6);
% %     mapData_orig = orig_Ez;
% %     vectorData_orig=orig_U(:,1:3); 
% %     coordinateOrigin = -1*X0(1:3);
% 
% %     C.Data = big_new_data_natural;
% %     C.Data = [data(:,1:3) big_new_data_natural(:,4:6)];
% %     orig_U = tdata(:,4:6)- repmat(C0(:,4:6),size(tdata,1),1);
%     d3D_col_orig = bsxfun(@minus, initial.Data,C0); %correcting RBD in initial data
% 
% %%%%%%%%%%
% OC5 - replaced by function: calcOrthoStrain
% 
% %% 5.1 Strain Calculation 3D
% %NEEDS TIDYING AND MAKING MORE GENERAL  for orthogonal directions
% 
% lx = dataSize(1);
% ly = dataSize(2);
% lz = dataSize(3);
% 
% tempUz1=zeros(lx,ly,lz+1);
% tempUz2=zeros(lx,ly,lz+1);
% tempUz1(:,:,1:lz)=orig_Uz;
% tempUz2(:,:,2:lz+1)=orig_Uz;
% % tempEz=-(tempUz1-tempUz2)./(nanmean(nanmean(orig_z(:,:,2)))-...
% %     nanmean(nanmean(orig_z(:,:,1))));
% tempEz=-(tempUz1-tempUz2)./dx_mm;
% 
% orig_Ez=tempEz(:,:,1:lz);
% 
% tempUz1=zeros(lx,ly,lz+1);
% tempUz2=zeros(lx,ly,lz+1);
% tempUz1(:,:,1:lz)=cleaned_Uz;
% tempUz2(:,:,2:lz+1)=cleaned_Uz;
% % tempEz=(tempUz1-tempUz2)./(nanmean(nanmean(cleaned_z(:,:,2)))-...
% %     nanmean(nanmean(cleaned_z(:,:,1))));
% tempEz=-(tempUz1-tempUz2)./dx_mm;
% 
% cleaned_Ez=tempEz(:,:,1:lz);
% %%%%%%%%%%%%%%%%%
% OC6 - where required individual functions reshape data rather than
% maintaining a global set
% 
% % This section is only required for rotation and strain - could be removed
% % cleaned_Ux=reshape(C.Data(:,4),C.Size);
% % cleaned_Uy=reshape(C.Data(:,5),C.Size);
% % cleaned_Uz=reshape(C.Data(:,6),C.Size);
% % 
% % orig_x=reshape(tdata(:,1),dataSize);
% % orig_y=reshape(tdata(:,2),dataSize);
% % orig_z=reshape(tdata(:,3),dataSize);
% % 
% % 
% % orig_Ux=reshape(orig_U(:,1),dataSize);
% % orig_Uy=reshape(orig_U(:,2),dataSize);
% % orig_Uz=reshape(orig_U(:,3),dataSize);
% %%%%%%%%%%%%%%%%
% OC7 - converted to xbinData
% 
% %     binSize = floor(C.Size(1)/numberOfBins);
% %     crUx_xbin = zeros(numberOfBins,C.Size(2),C.Size(3));
% %     crUy_xbin = crUx_xbin;
% %     crUz_xbin = crUx_xbin;
% % 
% %     for i = 1:numberOfBins
% %         fSlice = (i-1)*binSize+1;
% %         lSlice = i*binSize;
% %         crUx_xbin(i,:,:) = nanmean(cleaned_Ux_rot(fSlice:lSlice,:,:),1);
% %         crUy_xbin(i,:,:) = nanmean(cleaned_Uy_rot(fSlice:lSlice,:,:),1);
% %         crUz_xbin(i,:,:) = nanmean(cleaned_Uz_rot(fSlice:lSlice,:,:),1);
% %     end
% %     
% % % %     Ux2DAve = zeros(numberOfBins,1);
% % %     Ux2DAve = nanmean(nanmean(cleaned_Uy_rot(fSlice:lSlice,:,:)));
% % %     plot(squeeze(Ux2DAve)) 
% %     midX = round((C.Size(1:2)-1)/2)*dxOut;
% %     xint =2*midX(1)/numberOfBins;
% %     xHi = midX-0.5*xint;
% %     [ybin xbin zbin] = meshgrid(-midX(2):C.DX:midX(2),-xHi(1):xint:xHi(1),(1:lz)*C.DX);
% %     data3D_col = [xbin(:) ybin(:) zbin(:) crUx_xbin(:) crUy_xbin(:) crUz_xbin(:)];
% %%%%%%%%%%%%%%%
% OC8 - converted to function imRotateShell
% 
% if GeometryReference && mod(theta_deg,360) ~= 0
% 
% %     cleaned_Ux_rot = imrotateCustom(cleaned_Ux,theta_deg);
% %     cleaned_Uy_rot = imrotateCustom(cleaned_Uy,theta_deg);
% %     cleaned_Uz_rot = imrotateCustom(cleaned_Uz,theta_deg);
% % %     dataSize_rot = size(cleaned_Ux_rot);
% %     C.Size =size(cleaned_Ux_rot);
% % %     disp('WARNING: Check sign of theta_deg and trig in S4.3')
% %     disp(['Size of rotated data /voxel: (lx ly lz) = (' num2str(dataSize_rot) ').']); 
% %     disp(' ');
% %     
% %     %%
% %     if mod(theta_deg, 90) == 0
% % %         dx_rot = dx_mm;
% %     else
% %         C.DX = dx_mm/cosd(theta_deg);      %check trig!
% %     end
% %     
% %     midX = round((C.Size(1:2)-1)/2)*C.DX;
% %     [yrot xrot zrot] = meshgrid(-midX(2):dx_rot:midX(2),-midX(1):dx_rot:midX(1),(1:lz)*dx_rot);
%     %{
%     coordData = [xrot(:) yrot(:) zrot(:)];
%     mapData = [];
%     vectorData =[cleaned_Ux_rot(:) cleaned_Uy_rot(:) cleaned_Uz_rot(:)];
%     Title = {true;'(geo. test)'; cycleName};
% %     coordinateOrigin = [0 0 0];
% %     displacementsOrigin = [0 0 0];
%     for i=1:size(Plots,1)
%         mapID = Plots{i,1}; %[2 Y Map Vectors]; 
%         disp(['Plot ' num2str(i) ': {[' num2str(mapID) '] ' Plots{i,2} '}...']);
%         figure
%         plotOrthoSlice(mapID,coordData,dataSize_rot,mapData,vectorData,...
%             coordinateOrigin,displacementsOrigin,outOpts,Title);
%     end
%     %}
% else
% %     cleaned_Ux_rot = cleaned_Ux;
% %     cleaned_Uy_rot = cleaned_Uy;
% %     cleaned_Uz_rot = cleaned_Uz;
% % %     dataSize_rot = dataSize;
% % %     C.Size = dataSize_initial
% % %     dx_rot = dx_mm;
% %     midX = round((dataSize_rot(1:2)-1)/2)*dx_rot;
% end
% 
% %%%%%%%%%%
% OC9
%             if ExtractStrain && (mapID(3) <=6 || mapID(3) >=4)  
%                 mapID_bin = [mapID_bin(1:2) mapID(3)-3 0];
%                 plotTitle = {true;'(prebin)';cycleName};
%                 figure
%                 plotOrthoSlice(mapID_bin,binStrain.Data,binStrain.Size,...
%                     coordOrigin,displacementsOrigin,outOptsBin,plotTitle,binStrain.DX); 
%             end
% %}