function [ptPos,fields] = readVolvtk_UnStructGrid(varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Read the point-cloud (unstructured) vtk data 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [ptPos,fields] = readVolvtk_UnStructGrid(nfich)
%   [ptPos,fields] = readVolvtk_UnStructGrid(nfich,nfields)
% ----------------------------------------------------
%
%   Inputs:
%       - nfich: file name to be read
%       - nfields: number of fields expected to be read
%
%   Outputs:
%       - ptPos: the position coordiantes (x,y,z) of points in the vtk file
%       - fields: all fields in the vtk file (structure type)
% ================================================================
%
%   10/02/2016 Y.Chen
%

narginchk(1,2); % Check the number of input arguments.
nfich=varargin{1}; nfields=[];
disp(['reading the UnStructGrid vtk-file: ',nfich]);
disp '---------------'
if nargin==2
    nfields=varargin{2};
    disp 'You must be sure about the number of fields in the file to be read, otherwise, there will be no end !'
end;
% =======================================
fid=fopen(nfich,'r');
% Data format (binary/ASCII)
for i=1:3
    format=fgetl(fid);
end
disp(['reading ',format, ' data: ']);
numline=3;
% Point number
A='aaaaaa';
while (strcmp(A,'POINTS')==0)
    ligne = fgetl(fid);
    if (size(ligne,2)>=6);  A = ligne(1:6);  end
    numline=numline+1;
end
npts = sscanf(ligne(8:end),'%d');
typeVTK = sscanf(ligne(8+numel(num2str(npts))+1:end),'%s');
disp(['reading the position coordinates ...']);tic
% Read the point positions
if strcmp(format,'BINARY')
    typeMAT=typeVTK2MAT(typeVTK);
    B=fread(fid,npts*3,typeMAT,0,'b');
    ptPos = zeros(3,npts,typeMAT);
    ptPos(:)=B; ptPos=ptPos';
elseif strcmp(format,'ASCII')
    ptPos = dlmread(nfich,'',[numline 0 numline+npts-1 2]);
end; elpsT=toc;
    disp(['           Time used: ',num2str(elpsT),' s']);

% Read fields
if isempty(nfields); nfields=0; disp('only positions are read'); end
fields=cell(nfields,1); fieldnames=cell(nfields,1);

for i_field=1:nfields
    disp(['reading the ',num2str(i_field),'-th field: ']);
    A='aaaaaaa';
    while (strcmp(A,'SCALARS')==0)
        ligne = fgetl(fid);
        if (size(ligne,2)>=7);  A = ligne(1:7);  end
        numline=numline+1;
    end;
    fieldnames{i_field} = sscanf(ligne(9:end),'%c');
    typeVTK = fieldnames{i_field}(find(fieldnames{i_field}==' ')+1:end);
    fieldnames{i_field} = fieldnames{i_field}(1:find(fieldnames{i_field}==' ')-1);
%     disp(['    ',fieldnames{i_field},' ...']);
    A='aaaaaa';
    while (strcmp(A,'LOOKUP')==0)
        ligne = fgetl(fid);
        if (size(ligne,2)>=6);  A = ligne(1:6);  end
        numline=numline+1;
    end;
    if strcmp(format,'BINARY')
        typeMAT = typeVTK2MAT(typeVTK);
        fieldI=fread(fid,npts,typeMAT,0,'b');
    elseif strcmp(format,'ASCII')
        fieldI = single(dlmread(nfich,'',[numline 0 numline+npts-1 0]));
    end
    fields{i_field} = fieldI;
        elpsT=toc;
        disp([fieldnames{i_field},'(',typeVTK,typeMAT,')','     Time used: ',num2str(elpsT),' s']);
end

% regroup the fields into structure data
fields = cell2struct(fields,fieldnames,1);


fclose(fid);


function [typeMAT] = typeVTK2MAT(typeVTK)
if (strcmp(typeVTK,'bit')==1)
    typeMAT = 'ubit1';
elseif (strcmp(typeVTK,'unsigned_char')==1)
    typeMAT = 'uint8';
elseif (strcmp(typeVTK,'char')==1)
    typeMAT = 'int8';
elseif (strcmp(typeVTK,'unsigned_short')==1)
    typeMAT = 'uint16';
elseif (strcmp(typeVTK,'short')==1)
    typeMAT = 'int16';
elseif (strcmp(typeVTK,'unsigned_int')==1)
    typeMAT = 'uint32';
elseif (strcmp(typeVTK,'int')==1)
    typeMAT = 'int32';
elseif (strcmp(typeVTK,'unsigne_long')==1)
    typeMAT = 'uint64';
elseif (strcmp(typeVTK,'long')==1)
    typeMAT = 'int64';
elseif (strcmp(typeVTK,'float')==1)
    typeMAT = 'single';
elseif (strcmp(typeVTK,'double')==1)
    typeMAT = 'double';
end
end

end