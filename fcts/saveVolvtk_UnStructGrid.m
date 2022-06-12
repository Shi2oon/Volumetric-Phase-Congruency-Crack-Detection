function saveVolvtk_UnStructGrid(varargin)
% save a 3d point cloud into vtk format
%
%   saveVolvtk_UnStructGrid(fname,coords,name1,v1);
%   saveVolvtk_UnStructGrid(fname,coords,name1,v1,name2,v2,...);
% --------------------------------------------------------------------------------------------
%
%   Inputs:
%       - fname: file name
%       - coords: coordiantes of all input voxels (coords=[x,y,z])
%       - nameX: name of field
%       - vX: values of field
% =====================================================================
%  Attention: Not all data types are supported for all classes. (tested manually)
%  Here, the coordinates can be uint16, int16, uint32, int32, single and double
%        the component values can be uint8, int8, uint16, int16, uint32, int32, int64, single and double
% =======================================================================================================
%   Y.CHEN (10/2015)


% narginchk(3,14); % Check the number of input arguments.
fname = varargin{1};
disp(['saving the points cloud data in vtk-file: ',fname]);
disp 'ATENTION : Ensure that there is no NaN value in your inputs !'
tic
% Header of the vtk file for unstructured grid data
% -------------------------------------------------
fid = fopen(fname, 'w');
fprintf(fid, '# vtk DataFile Version 4.5\n');
fprintf(fid, 'tomo_Volume\n');
fprintf(fid, 'BINARY\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

% Point positions (x y z)
% -----------------------
coords = varargin{2}; npts = size(coords,1);
typeMAT=class(coords);  [typeVTK,typeMAT]=type_MAT2VTK(typeMAT);
fprintf(fid, ['POINTS ',num2str(npts),' ',typeVTK,'\n']);
fwrite(fid,coords',typeMAT,0,'b'); % fwrite is in colomn order, while paraview read in row order
disp(['Position data: ',typeMAT,typeVTK]);

% Cell indice and Cell type
% -------------------------------
ncells = npts;   % each point is one cell
fprintf(fid, ['\nCELLS ',num2str(ncells),' ',num2str(npts+ncells),' \n']);
fwrite(fid,[ones(1,ncells);[1:npts]],'int32',0,'b');   % the cell types and indices must be int format

fprintf(fid, ['\nCELL_TYPES ',num2str(ncells),' \n']);
fwrite(fid,ones(ncells,1),'int32',0,'b');   % the cell types and indices must be int format
    
    i0=2; nfields = (numel(varargin)-i0)/2; % number of fields to be saved


% Field data
% -----------
fprintf(fid, ['\nPOINT_DATA ',num2str(npts),' \n']);
if nfields-double(uint64(nfields))~=0;  disp('missing input information !');fclose(fid);return; end
for ifield=1:nfields
    nameI = varargin{i0+(ifield-1)*2+1};
    fieldI = varargin{i0+ifield*2};
    typeMAT=class(fieldI);   [typeVTK,typeMAT]=type_MAT2VTK(typeMAT);
    fprintf(fid, ['SCALARS ', nameI,' ',typeVTK,'\n']);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fwrite(fid,fieldI(:),typeMAT,0,'b');
    disp(['Field data ',nameI,': ',typeMAT,typeVTK])
end

fclose(fid);
toc

end

function [typeVTK,typeMAT] =  type_MAT2VTK(typeMAT)
if (strcmp(typeMAT,'logical')==1)
    typeVTK = 'bit'; typeMAT = 'ubit1';
elseif (strcmp(typeMAT,'uint8')==1)
    typeVTK = 'unsigned_char';
elseif (strcmp(typeMAT,'int8')==1)
    typeVTK = 'char';
elseif (strcmp(typeMAT,'uint16')==1)
    typeVTK = 'unsigned_short';
elseif (strcmp(typeMAT,'int16')==1)
    typeVTK = 'short';
elseif (strcmp(typeMAT,'uint32')==1)
    typeVTK = 'unsigned_int';
elseif (strcmp(typeMAT,'int32')==1)
    typeVTK = 'int';
elseif (strcmp(typeMAT,'uint64')==1)
    typeVTK = 'unsigne_long';
elseif (strcmp(typeMAT,'int64')==1)
    typeVTK = 'long';
elseif (strcmp(typeMAT,'single')==1)
    typeVTK = 'float';  typeMAT = 'float';
elseif (strcmp(typeMAT,'double')==1)
    typeVTK = 'double';
end

end


