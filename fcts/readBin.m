function [X] = readBin(varargin)
% read bin file 
%   [X] = readBin(fname)
%   [X] = readBin(fname,nblocs)
%   ---------------------------
%
%   Inputs
%       >  fname : file name to be read
%       > nblocs : (optional) define the number of blocs for reading
%
%   Outputs
%       >      X : data in fname
%
% Modified by Yang CHEN 2018.02.24
%          add an option to read multicolomn matrix
%

disp(['reading bin-file: ',varargin{1}])

fid = fopen(varargin{1});

% reading blocs
if nargin==2
    nblocs = varargin{2};
else
    nblocs = 1;
end

% points number
npts = fgetl(fid);
npts = str2num(npts);
sizX = [];
if length(npts)>1
    sizX = npts;
    npts = prod(npts);
end

% data type
typeVTK = strtrim(fgetl(fid));

typeMAT = typeVTK2MAT(typeVTK);

% reading
inc=ceil(npts/nblocs); over=0;
if strcmp(typeMAT,'ubit1')
    X = false(npts,1);
else
    X = zeros(npts,1,typeMAT);
end
for ibloc = 1:nblocs
    i1 = 1+inc*(ibloc-1);
    i2 = i1+inc-1;
    if i2>npts; i2=npts; inc=i2-i1+1; over=1;  end
    X(i1:i2)=fread(fid,inc,typeMAT,0,'b');

    disp(['bloc',num2str(ibloc),': ',num2str([i1,i2,inc])])
    if over==1; break; end
end

fclose(fid);

% if the data is multicolomn matrix
if ~isempty(sizX)
    X = reshape(X,sizX);
end

end

%% =======================================
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
elseif (strcmp(typeVTK,'unsigned_long')==1)
    typeMAT = 'uint64';
elseif (strcmp(typeVTK,'long')==1)
    typeMAT = 'int64';
elseif (strcmp(typeVTK,'float')==1)
    typeMAT = 'single';
elseif (strcmp(typeVTK,'double')==1)
    typeMAT = 'double';
end
end