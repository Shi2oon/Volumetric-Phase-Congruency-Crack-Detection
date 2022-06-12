function []=LTrack_fromStart(CC1,CC2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Track cracks from step1 and give them the same labels for the following
%%% steps  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [connectL_sort_1in2]=LsortTracking_btw2steps(fname1,fname2,sizVOL);
%   [connectL_sort_1in2]=LsortTracking_btw2steps(...,option);
%   ============================================================
% Inputs:
%   > fname1,fname2: the file names of the two crack volumes
%               Note, the format of input crack volumes need to respect a
%               special format (see saveVolvtk_UnStructGrid.m)
%   > sizVOL: the size of the volume to be re-built
%   > option: 'remain' or 'all' indicating whether only remaining cracks or
%               both remaining and newborn cracks are tracked
%   -----------------------------------------------------------------------
% Outputs:
%   > Ltrack: the labels of (1) remaining cracks in step2 with the sorting
%   number of the step1, (2) new-born cracks in step2
%   ptPos: the positions of points in step2
%   fields2: the field information in step2
%   -----------------------------------------------------------------------
%
% 14/03/2016 Y.Chen
%


sizVOL=varargin{3};
option=varargin{4};

tic
nR=sizVOL(1);  nC=sizVOL(2);  nB=sizVOL(3);

% ----- read the first crack volume
[ptPos1,fields1]=readVolvtk_UnStructGrid(fname1,5);
VxLst1=sub2ind([nR,nC,nB],ptPos1(:,2),ptPos1(:,1),ptPos1(:,3));

% ----- read the second crack volume
[ptPos2,fields2]=readVolvtk_UnStructGrid(fname2,4);
VxLst2=sub2ind([nR,nC,nB],ptPos2(:,2),ptPos2(:,1),ptPos2(:,3));
npts_2=length(VxLst2);

% ----- tracking the remaining cracks
Lid_remain_2 = fields2.connectL_sort(ismember(VxLst2,VxLst1));
Lid_remain_2_uniq = unique(Lid_remain_2);
% label the remaining cracks with the label ids of the previous step
if strcmp(option,'remain')
    dataType = GetType(max(fields1.Ltrack));
elseif strcmp(option,'all')
    dataType = GetType(max(fields1.Ltrack)+max(fields2.connectL_sort));
end
Ltrack = zeros(npts_2,1,dataType);
for i=1:numel(Lid_remain_2_uniq)
    idx_remain_2_i = fields2.connectL_sort == Lid_remain_2_uniq(i);
    idx_remain_1_i=ismember(VxLst1,VxLst2(idx_remain_2_i));
    Ltrack(idx_remain_2_i)=mode(fields1.Ltrack(idx_remain_1_i));                    
end
elasT=toc; disp(['Remaining labels tracking is completed. Time used: ',num2str(elasT)]);

if strcmp(option,'remain');
    disp('The new-born cracks are not labeled');return;
end


% ----- label the new-born cracks 
nbL2=length(unique(fields2.connectL_sort));
Lid_born_2 = setdiff([1:nbL2],Lid_remain_2_uniq);
Lid_endOF1 = max(fields1.Ltrack(ismember(VxLst1,VxLst2)));

for i=1:numel(Lid_born_2)  % give the label id to each new-born group in def2
    idx_born_2_i = fields2.connectL_sort == Lid_born_2(i);
    Ltrack(idx_born_2_i) = Lid_endOF1+i;
end

Ltrack = typeReduc(Ltrack);

elasT=toc; disp(['New-born labels tracking is completed. Time used: ',num2str(elasT)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataType]=GetType(valueMax)
    if valueMax <= intmax('uint8')
        dataType = 'uint8';
    elseif valueMax <= intmax('uint16')
        dataType = 'uint16';
    elseif valueMax <= intmax('uint32')
        dataType = 'uint32';
    else
        dataType = 'double';
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X] = typeReduc(X)
    type=GetType(max(X(:)));
    if strcmp(type,'uint8')
        X=uint8(X);
    elseif strcmp(type,'uint16')
        X=uint16(X);
    elseif strcmp(type,'uint32')
        X=uint32(X);
    end  
