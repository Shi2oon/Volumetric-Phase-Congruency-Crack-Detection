function [Wav,WavS,S,L,W,H,G,V,dS,Vav,cONa,bONa] = ...
CrackMEASURE_2(VxLst,VOLres,Vgap,config,lambda,nsub,Ot,ex,ey,ez,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%           measurement of the average opening of each crack
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Y.Chen 20/07/2016

if nargin==10
    method = 'half';
elseif nargin==11
    method = varargin{1};
end


sizVOL= size(VOLres);

%% compute the air fraction within each crack voxel
VxLst_array=cell2mat(VxLst);
Vres = VOLres(VxLst_array);
Qair = Qair_comput(Vres,Vgap(VxLst_array),intmax(class(VOLres)),method);

figure;hist(Qair(:),1000)

Qair(Qair>1)=1; Qair(Qair<0)=0;

% define a volume describing the air fraction at each crack voxel
VOL_Qair = zeros(sizVOL,'single');
VOL_Qair(VxLst_array) = Qair;


%% Measurements
if iscell(VxLst)
    ncracks = length(VxLst);

    % allocate the memory
    Wav = zeros(ncracks,1);
    WavS = zeros(ncracks,1);
    S = zeros(ncracks,1);
    dS = cell(ncracks,1);
    Vav = cell(ncracks,1);
    
    L = zeros(ncracks,1);
    W = zeros(ncracks,1);
    H = zeros(ncracks,1);
    G = zeros(ncracks,3);
    V = zeros(ncracks,1);

    cONa = zeros(ncracks,1);
    bONa = zeros(ncracks,1);
 
    %					              loop for each crack
    % --------------------------------------------------------------------
    for icrack = 1:ncracks
        VxLstI = VxLst{icrack};
       
        % ---- crack opening, surface area
        QairI = VOL_Qair(VxLstI);
        VOL_QairI = zeros(sizVOL,'single');
        VOL_QairI(VxLstI) = QairI;
        %[Wav(icrack),S(icrack)] = surfThickness(VxLstI,VOL_QairI,...
        %                     config,lambda,nsub,style,binCapacity,minBin);
                     
        [Wav(icrack),WavS(icrack),S(icrack),dS{icrack},Vav{icrack}] = ...
                      surfThickness_2(VxLstI,VOL_QairI,config,lambda,nsub);
                                                         
        % ---- crack dimension, location
        [L(icrack),W(icrack),H(icrack),G(icrack,:),V(icrack)] = ...
                                      objDimens(VxLstI,sizVOL,Ot,ex,ey,ez);

        % ---- compute aspect ratios of each crack in cylin coordinates
	    [cONa(icrack),bONa(icrack)] = ...
                            AspRatio_cylin(VxLstI,sizVOL,Ot,ex,ey,ez);
                
    end

elseif ~iscell(VxLst)
    % crack opening, surface area
	[Wav,S] = ...
	surfThickness_2(VxLst,VOL_Qair,config,lambda,nsub);

    % crack dimension, location
    [L,W,H,G,V] = objDimens0(VxLst,sizVOL,Ot,ex,ey,ez);

    % compute aspect ratios of each crack in cylin coordinates
    [cONa,bONa] = AspRatio_cylin(VxLst,sizVOL,Ot,ex,ey,ez);
 
end








